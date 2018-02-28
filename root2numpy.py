
import ROOT
from ROOT import AddressOf
import numpy as np
import time, datetime, sys, os, re

import scipy.weave
from timer import timer

cast_map = {
    'F': 'Float_t',
    #'I': 'ULong64_t',
    'I': 'Int_t',
    'B': 'Char_t',
    }
            

class root2numpy:
    
    def python( self, fname, keys = None ):
        t = timer("python")
        print 'parseTree'
        if not os.path.exists( fname ):
            print "catalog not found", fname
            return None
        if not '.root' in fname:
            print "file is not a root file", fname
            return None
        fstream = ROOT.TFile( fname )
        result = {}
        for key in fstream.GetListOfKeys():
            tree = fstream.Get( key.GetName() )
            print key.GetName()
            list_of_branches = [ branch.GetTitle().split('/') for branch in tree.GetListOfBranches() if len(branch.GetTitle().split('/')) == 2 ]# and len(branch.GetTitle().split('[')) == 1 ]
            #declaration_list = [ "%s %s;"%( cast_map[branch[1]], branch[0] ) for branch in list_of_branches ]
            declaration_list = [ "%s %s;"%( cast_map[branch[1]], re.sub("(.*?)\[([^\d]\w*?)\]", "\\1[100000]", branch[0]) ) for branch in list_of_branches ]
            branch_titles = zip(*list_of_branches)[0]
            struct_string = "struct struct_%s{%s};"%( key.GetName(),' '.join( declaration_list ))
            #print struct_string
            branch_scalar_names = [ title for title in branch_titles if len(title.split('[')) == 1 ]
            branch_array_names = [ title.split('[')[0] for title in branch_titles if len(title.split('[')) > 1 ]
            branch_array_names_sizes = [ [ title.split('[')[0], re.search('\[(.*?)\]', title).group(1) ] for title in branch_titles if len(title.split('[')) > 1 ]
            #for name, size in branch_array_names_sizes:
                #print name, size
            dictionary = {}
            ROOT.gROOT.ProcessLine( struct_string )
            entry = eval( "ROOT.struct_%s()"%key.GetName() )
            
            for name in branch_scalar_names + branch_array_names:
                tree.SetBranchAddress( name, ROOT.AddressOf( entry, name ) )
                dictionary[name]= []
            for index in xrange(tree.GetEntries()):
                tree.GetEntry(index)
                for name in branch_scalar_names:
                    dictionary[name].append( getattr( entry, name ) )
                    #dictionary[name].append( entry.__dict__[name] )
                for name, size in branch_array_names_sizes:
                    actual_size = (int(size) if size.isdigit() else getattr(entry, size) )
                    #dictionary[name].append( getattr( entry, name )[:actual_size] )
                    dictionary[name].append( [ getattr( entry, name )[pix] for pix in xrange( actual_size ) ] )
                    #dictionary[name].append( [ entry.__dict__[name][pix] for pix in xrange( actual_size ) ] )
                    #print name, getattr( entry, name )[:(int(size) if size.isdigit() else getattr(entry, size) )]
            
            tree.ResetBranchAddress(0)
            #for name in branch_array_names:
                #dictionary[name] = np.array( dictionary[name] )
            result[key.GetName()] = dictionary

        del t
        return result
    
    def weave( self, fname, keys = None ):
        t = timer("weave")
        print 'parseTree'
        if not os.path.exists( fname ):
            print "catalog not found", fname
            return None
        if not '.root' in fname:
            print "file is not a root file", fname
            return None
        fstream = ROOT.TFile( fname )
        result = {}
        for key in fstream.GetListOfKeys():
            tree = fstream.Get( key.GetName() )
            print key.GetName()
            list_of_branches = [ branch.GetTitle().split('/') for branch in tree.GetListOfBranches() if len(branch.GetTitle().split('/')) == 2 ]# and len(branch.GetTitle().split('[')) == 1 ]
            declaration_list = [ "%s %s;"%( cast_map[branch[1]], re.sub("(.*?)\[([^\d]\w*?)\]", "\\1[100000]", branch[0]) ) for branch in list_of_branches ]
            branch_titles = zip(*list_of_branches)[0]
            struct_string = "struct struct_%s{%s};"%( key.GetName(),' '.join( declaration_list ))
            #print struct_string
            branch_scalar_names = [ title for title in branch_titles if len(title.split('[')) == 1 ]
            branch_array_names = [ title.split('[')[0] for title in branch_titles if len(title.split('[')) > 1 ]
            branch_array_names_sizes = [ [ title.split('[')[0], re.search('\[(.*?)\]', title).group(1) ] for title in branch_titles if len(title.split('[')) > 1 ]
            dictionary = {}
            ROOT.gROOT.ProcessLine( struct_string )
            entry = eval( "ROOT.struct_%s()"%key.GetName() )
            
            for name in branch_scalar_names + branch_array_names:
                tree.SetBranchAddress( name, ROOT.AddressOf( entry, name ) )
                dictionary[name]= []

            for index in xrange(tree.GetEntries()):
                tree.GetEntry(index)
                for name in branch_scalar_names:
                    dictionary[name].append( getattr( entry, name ) )
                for name, size in branch_array_names_sizes:
                    actual_size = (int(size) if size.isdigit() else getattr(entry, size) )
                    code = """
                        py::list l(actual_size);
                        for(int i=0; i<actual_size; ++i){
                            l[i] = val[i];
                        }
                        current.append(l);
                        return_val = l;
                    """
                    val = getattr( entry, name )
                    current = dictionary[name]
                    scipy.weave.inline( code, ['actual_size', 'val', 'current'] )
                    #current.append( scipy.weave.inline( code, ['actual_size', 'val'] ) )
                    #dictionary[name].append( [ getattr( entry, name )[pix] for pix in xrange( actual_size ) ] )
            
            tree.ResetBranchAddress(0)
            #for name in branch_array_names:
                #dictionary[name] = np.array( dictionary[name] )
            result[key.GetName()] = dictionary

        t()
        return result


    def weave2( self, fname, keys = None ):
        t = timer("weave2")
        print 'parseTree'
        if not os.path.exists( fname ):
            print "catalog not found", fname
            return None
        if not '.root' in fname:
            print "file is not a root file", fname
            return None
        #code = """
            #std::cout << "name of file " << fname.c_str() << std::endl;
            #TFile f(fname.c_str());
            #TList *keys = f.GetListOfKeys();
            #TIter next_key(keys);
            #TObject *key = 0;
            #while( ( key = next_key() ) ){
                #std::cout << key->GetName() << std::endl;
                #TTree *t = (TTree*)f.Get( key->GetName() );
                #TObjArray *branches = t->GetListOfBranches();
                #std::cout << "list of branches " << branches->GetEntries() << std::endl;
                #std::string 
                #for( int i = 0; i < branches->GetEntries(); ++i){
                    #std::cout << (*branches)[i]->GetName() << std::endl;
                    #(*branches)[i]->Print();
                #}
                #int n = (int) t->GetEntries();
                #std::cout << "size of tree " << n << std::endl;
            #}
        #""" #% {
            ##'declarations' : "\n".join( [ "py::list %s;" %(name) for name in branch_scalar_names ] ),
            ##'assignments' : "\n".join( [ "%s.append( entry.%s );" %(name, name) for name in branch_scalar_names ] )
            ##}
        #print code
        #scipy.weave.inline( code, [ 'fname' ], 
                           #headers=['"TFile.h"', '"TTree.h"', '"TObject.h"', '"TObjArray.h"'], 
                           #libraries=['Core'],
                           #include_dirs=['/usr/include/root/'],
                           #library_dirs=['/usr/lib/x86_64-linux-gnu/root5.34/']
                           #)
        
        fstream = ROOT.TFile( fname )
        headers = {}
        result = {}
        for key in fstream.GetListOfKeys():
            tree = fstream.Get( key.GetName() )
            print key.GetName()
            unique_branches = set([ branch.GetName() for branch in tree.GetListOfBranches() ])
            #print unique_branches
            branches = [ tree.GetBranch(branch_name) for branch_name in unique_branches ]
            headers[key.GetName()] = [ {
                "type": cast_map[ branch.GetTitle().split('/')[-1] ] if len(branch.GetTitle().split('/')) > 1 else None, 
                 "name": branch.GetName(), 
                 "size": re.search('\[(.*?)\]', branch.GetTitle()).group(1) if '[' in branch.GetTitle() else None
                     } for branch in branches ]
            print branches
        fstream.Close()
        result = {}
        for key, branches in headers.items():
            if key == 'config': continue
            #print key
            size = len([ branch['name'] for branch in branches if branch['type'] ] )
            #print size
            declarations = '\n'.join( set([ 'py::list list_%(name)s; %(type)s %(name)s'%branch + ('[%s];'%( branch['size'] if branch['size'].isdigit() else '100000' ) if branch['size'] else ';')  for branch in branches if branch['type'] ]) )
            #print declarations
            addresses = '\n'.join( set([ 'tree->SetBranchAddress("%(name)s", &%(name)s );'%branch if not branch['size'] else 'tree->SetBranchAddress("%(name)s", %(name)s );'%branch for branch in branches if branch['type'] ]) )
            #print addresses
            assignments = '\n'.join( set([ 
                'list_%(name)s.append( %(name)s );'%branch 
                for branch in branches if branch['type'] and not branch['size'] 
                ]) )
            copies = '\n'.join( set([ 'py::tuple entry_%(name)s(2); entry_%(name)s[0] = "%(name)s"; entry_%(name)s[1] = list_%(name)s; '% branch +  'results[%s] = entry_%s;' % (index,branch['name']) for index, branch in zip(range(size),branches) if branch['type'] and not branch['size'] ]) )
            #print assignments
            code = """
std::cout << "name of file " << fname.c_str() << std::endl;
TFile f(fname.c_str());
TTree *tree = (TTree*) f.Get(key.c_str());
%(declarations)s
%(addresses)s
int n = (int) tree->GetEntries();
std::cout << "entries " << n  << std::endl;
for(int i=0;i<n;++i){
tree->GetEntry(i);
%(assignments)s
}
py::tuple results(%(size)s);
%(copies)s
return_val = results;

            """ %{ 
                'declarations': declarations, 
                'addresses': addresses,
                'assignments': assignments,
                'size': size,
                'copies': copies,
                }
            print code
            partial_result = scipy.weave.inline( code, [ 'fname', 'key' ], 
                           headers=['"TFile.h"', '"TTree.h"', '"TObject.h"', '"TObjArray.h"'], 
                           libraries=['Core'],
                           include_dirs=['/usr/include/root/'],
                           library_dirs=['/usr/lib/x86_64-linux-gnu/root5.34/'],
                           #compiler= 'gcc -O3'
                           )
            #print partial_result
            #result[key] = dict(partial_result)
        t()
        return result

    def weaveFull( self, fname ):
        t = timer("weave2")
        print 'parseTree'
        if not os.path.exists( fname ):
            print "catalog not found", fname
            return None
        if not '.root' in fname:
            print "file is not a root file", fname
            return None

        code = """
std::cout << "name of file " << fname.c_str() << std::endl;
TFile f(fname.c_str());
TTree *tree = (TTree*) f.Get("hitSumm");
int n = (int) tree->GetEntries();
py::tuple list_flag(n); Int_t flag;
py::tuple list_ohdu(n); Int_t ohdu;
py::tuple list_runID(n); Int_t runID;
py::tuple list_E0(n); Float_t E0;
py::tuple list_n0(n); Float_t n0;
py::tuple list_nSavedPix(n); Int_t nSavedPix;
py::tuple list_xPix(n); Int_t xPix[100000];
py::tuple list_yPix(n); Int_t yPix[100000];
py::tuple list_ePix(n); Float_t ePix[100000];
py::tuple list_level(n); Int_t level[100000];
tree->SetBranchAddress("flag", &flag );
tree->SetBranchAddress("ohdu", &ohdu );
tree->SetBranchAddress("runID", &runID );
tree->SetBranchAddress("E0", &E0 );
tree->SetBranchAddress("n0", &n0 );
tree->SetBranchAddress("nSavedPix", &nSavedPix );
tree->SetBranchAddress("xPix", xPix );
tree->SetBranchAddress("yPix", yPix );
tree->SetBranchAddress("level", level );
tree->SetBranchAddress("ePix", ePix );
std::cout << "entries " << n  << std::endl;
for(int i=0;i<n;++i){
tree->GetEntry(i);
list_flag[i] = flag;
list_ohdu[i] = ohdu;
list_runID[i] = runID;
list_E0[i] = E0;
list_n0[i] = n0;
list_nSavedPix[i] = nSavedPix;
py::tuple xPix_temp(nSavedPix);
py::tuple yPix_temp(nSavedPix);
py::tuple ePix_temp(nSavedPix);
py::tuple level_temp(nSavedPix);
for(int n=0;n<nSavedPix;++n){
xPix_temp[n] = xPix[n];
yPix_temp[n] = yPix[n];
ePix_temp[n] = ePix[n];
level_temp[n] = level[n];
}
list_xPix[i] = xPix_temp;
list_yPix[i] = yPix_temp;
list_ePix[i] = ePix_temp;
list_level[i] = level_temp;
}
py::tuple results(10);
py::tuple entry_flag(2); entry_flag[0] = "flag"; entry_flag[1] = list_flag; 
results[0] = entry_flag;
py::tuple entry_ohdu(2); entry_ohdu[0] = "ohdu"; entry_ohdu[1] = list_ohdu; 
results[1] = entry_ohdu;
py::tuple entry_runID(2); entry_runID[0] = "runID"; entry_runID[1] = list_runID; 
results[2] = entry_runID;
py::tuple entry_E0(2); entry_E0[0] = "E0"; entry_E0[1] = list_E0; 
results[3] = entry_E0;
py::tuple entry_n0(2); entry_n0[0] = "n0"; entry_n0[1] = list_n0; 
results[4] = entry_n0;
py::tuple entry_nSavedPix(2); entry_nSavedPix[0] = "nSavedPix"; entry_nSavedPix[1] = list_nSavedPix;
results[5] = entry_nSavedPix;
py::tuple entry_xPix(2); entry_xPix[0] = "xPix"; entry_xPix[1] = list_xPix;
results[6] = entry_xPix;
py::tuple entry_yPix(2); entry_yPix[0] = "yPix"; entry_yPix[1] = list_yPix;
results[7] = entry_yPix;
py::tuple entry_ePix(2); entry_ePix[0] = "ePix"; entry_ePix[1] = list_ePix;
results[8] = entry_ePix;
py::tuple entry_level(2); entry_level[0] = "level"; entry_level[1] = list_level;
results[9] = entry_level;

return_val = results;
            """
        #print code
        partial_result = scipy.weave.inline( code, [ 'fname' ], 
                        headers=['"TFile.h"', '"TTree.h"', '"TObject.h"', '"TObjArray.h"'], 
                        libraries=['Core'],
                        include_dirs=['/usr/include/root/'],
                        library_dirs=['/usr/lib/x86_64-linux-gnu/root5.34/'],
                        #compiler= 'gcc -O3'
                        extra_compile_args=['-O3'],
                        )
        del t
        return dict(partial_result)

    def weaveFullExtra( self, fname, selection = '' ):
        t = timer("weaveFullExtra")
        print 'parseTree'
        if not os.path.exists( fname ):
            print "catalog not found", fname
            return None
        if not '.root' in fname:
            print "file is not a root file", fname
            return None
        if not 'Extra.root' in fname:
            print "file is not a Extra root file", fname
            return None

        code = """
std::cout << "name of file " << fname.c_str() << std::endl;
TFile f(fname.c_str());
TTree *tree = (TTree*) f.Get("hitSumm");
int n = (int) tree->GetEntries();
py::tuple list_flag(n); Int_t flag;
py::tuple list_ohdu(n); Int_t ohdu;
py::tuple list_runID(n); Int_t runID;
py::tuple list_E0(n); Float_t E0;
py::tuple list_n0(n); Float_t n0;
py::tuple list_nSavedPix(n); Int_t nSavedPix;
py::tuple list_xPix(n); Int_t xPix[100000];
py::tuple list_yPix(n); Int_t yPix[100000];
py::tuple list_ePix(n); Float_t ePix[100000];
py::tuple list_level(n); Int_t level[100000];
py::tuple list_T(n); Float_t T;
py::tuple list_L(n); Float_t L;
py::tuple list_stdT(n); Float_t stdT;
py::tuple list_stdL(n); Float_t stdL;
py::tuple list_errT(n); Float_t errT;
py::tuple list_errL(n); Float_t errL;
tree->SetBranchAddress("flag", &flag );
tree->SetBranchAddress("ohdu", &ohdu );
tree->SetBranchAddress("runID", &runID );
tree->SetBranchAddress("E0", &E0 );
tree->SetBranchAddress("n0", &n0 );
tree->SetBranchAddress("nSavedPix", &nSavedPix );
tree->SetBranchAddress("xPix", xPix );
tree->SetBranchAddress("yPix", yPix );
tree->SetBranchAddress("level", level );
tree->SetBranchAddress("ePix", ePix );
tree->SetBranchAddress("T", &T );
tree->SetBranchAddress("L", &L );
tree->SetBranchAddress("stdT", &stdT );
tree->SetBranchAddress("stdL", &stdL );
tree->SetBranchAddress("errT", &errT );
tree->SetBranchAddress("errL", &errL );
std::cout << "entries " << n  << std::endl;
for(int i=0;i<n;++i){
tree->GetEntry(i);
list_flag[i] = flag;
list_ohdu[i] = ohdu;
list_runID[i] = runID;
list_E0[i] = E0;
list_n0[i] = n0;
list_nSavedPix[i] = nSavedPix;

py::tuple xPix_temp(nSavedPix);
py::tuple yPix_temp(nSavedPix);
py::tuple ePix_temp(nSavedPix);
py::tuple level_temp(nSavedPix);
for(int n=0;n<nSavedPix;++n){
xPix_temp[n] = xPix[n];
yPix_temp[n] = yPix[n];
ePix_temp[n] = ePix[n];
level_temp[n] = level[n];
}
list_xPix[i] = xPix_temp;
list_yPix[i] = yPix_temp;
list_ePix[i] = ePix_temp;
list_level[i] = level_temp;
list_T[i] = T;
list_L[i] = L;
list_stdT[i] = stdT;
list_stdL[i] = stdL;
list_errT[i] = errT;
list_errL[i] = errL;
}
py::tuple results(16);
py::tuple entry_flag(2); entry_flag[0] = "flag"; entry_flag[1] = list_flag; 
results[0] = entry_flag;
py::tuple entry_ohdu(2); entry_ohdu[0] = "ohdu"; entry_ohdu[1] = list_ohdu; 
results[1] = entry_ohdu;
py::tuple entry_runID(2); entry_runID[0] = "runID"; entry_runID[1] = list_runID; 
results[2] = entry_runID;
py::tuple entry_E0(2); entry_E0[0] = "E0"; entry_E0[1] = list_E0; 
results[3] = entry_E0;
py::tuple entry_n0(2); entry_n0[0] = "n0"; entry_n0[1] = list_n0; 
results[4] = entry_n0;
py::tuple entry_nSavedPix(2); entry_nSavedPix[0] = "nSavedPix"; entry_nSavedPix[1] = list_nSavedPix;
results[5] = entry_nSavedPix;
py::tuple entry_xPix(2); entry_xPix[0] = "xPix"; entry_xPix[1] = list_xPix;
results[6] = entry_xPix;
py::tuple entry_yPix(2); entry_yPix[0] = "yPix"; entry_yPix[1] = list_yPix;
results[7] = entry_yPix;
py::tuple entry_ePix(2); entry_ePix[0] = "ePix"; entry_ePix[1] = list_ePix;
results[8] = entry_ePix;
py::tuple entry_level(2); entry_level[0] = "level"; entry_level[1] = list_level;
results[9] = entry_level;
py::tuple entry_T(2); entry_T[0] = "T"; entry_T[1] = list_T;
results[10] = entry_T;
py::tuple entry_L(2); entry_L[0] = "L"; entry_L[1] = list_L;
results[11] = entry_L;
py::tuple entry_stdT(2); entry_stdT[0] = "stdT"; entry_stdT[1] = list_stdT;
results[12] = entry_stdT;
py::tuple entry_stdL(2); entry_stdL[0] = "stdL"; entry_stdL[1] = list_stdL;
results[13] = entry_stdL;
py::tuple entry_errT(2); entry_errT[0] = "errT"; entry_errT[1] = list_errT;
results[14] = entry_errT;
py::tuple entry_errL(2); entry_errL[0] = "errL"; entry_errL[1] = list_errL;
results[15] = entry_errL;

return_val = results;
            """
        #print code
        partial_result = scipy.weave.inline( code, [ 'fname' ], 
                        headers=['"TFile.h"', '"TTree.h"', '"TObject.h"', '"TObjArray.h"'], 
                        libraries=['Core'],
                        include_dirs=['/usr/include/root/'],
                        library_dirs=['/usr/lib/x86_64-linux-gnu/root5.34/'],
                        #compiler= 'gcc -O3'
                        extra_compile_args=['-O3'],
                        )
        del t
        return dict(partial_result)

    def weaveFullExtraSelection( self, fname, selection = '' ):
        t = timer("weaveFullExtraSelection")
        print 'parseTree'
        if not os.path.exists( fname ):
            print "catalog not found", fname
            return None
        if not '.root' in fname:
            print "file is not a root file", fname
            return None
        if not 'Extra.root' in fname:
            print "file is not a Extra root file", fname
            return None

        code = """
std::cout << "name of file " << fname.c_str() << std::endl;
TFile f(fname.c_str());
TTree *tree = (TTree*) f.Get("hitSumm");
int n = (int) tree->GetEntries();
py::list list_flag; Int_t flag;
py::list list_ohdu; Int_t ohdu;
py::list list_runID; Int_t runID;
py::list list_E0; Float_t E0;
py::list list_n0; Float_t n0;
py::list list_nSavedPix; Int_t nSavedPix;
py::list list_xPix; Int_t xPix[100000];
py::list list_yPix; Int_t yPix[100000];
py::list list_ePix; Float_t ePix[100000];
py::list list_level; Int_t level[100000];
py::list list_T; Float_t T;
py::list list_L; Float_t L;
py::list list_stdT; Float_t stdT;
py::list list_stdL; Float_t stdL;
py::list list_errT; Float_t errT;
py::list list_errL; Float_t errL;
tree->SetBranchAddress("flag", &flag );
tree->SetBranchAddress("ohdu", &ohdu );
tree->SetBranchAddress("runID", &runID );
tree->SetBranchAddress("E0", &E0 );
tree->SetBranchAddress("n0", &n0 );
tree->SetBranchAddress("nSavedPix", &nSavedPix );
tree->SetBranchAddress("xPix", xPix );
tree->SetBranchAddress("yPix", yPix );
tree->SetBranchAddress("level", level );
tree->SetBranchAddress("ePix", ePix );
tree->SetBranchAddress("T", &T );
tree->SetBranchAddress("L", &L );
tree->SetBranchAddress("stdT", &stdT );
tree->SetBranchAddress("stdL", &stdL );
tree->SetBranchAddress("errT", &errT );
tree->SetBranchAddress("errL", &errL );
std::cout << "entries " << n  << std::endl;
for(int i=0;i<n;++i){
tree->GetEntry(i);
%s
list_flag.append( flag );
list_ohdu.append( ohdu );
list_runID.append( runID );
list_E0.append( E0 );
list_n0.append( n0 );
list_nSavedPix.append( nSavedPix );
py::tuple xPix_temp(nSavedPix);
py::tuple yPix_temp(nSavedPix);
py::tuple ePix_temp(nSavedPix);
py::tuple level_temp(nSavedPix);
for(int n=0;n<nSavedPix;++n){
xPix_temp[n] = xPix[n];
yPix_temp[n] = yPix[n];
ePix_temp[n] = ePix[n];
level_temp[n] = level[n];
}
list_xPix.append( xPix_temp );
list_yPix.append( yPix_temp );
list_ePix.append( ePix_temp );
list_level.append( level_temp );
list_T.append( T );
list_L.append( L );
list_stdT.append( stdT );
list_stdL.append( stdL );
list_errT.append( errT );
list_errL.append( errL );
}
py::tuple results(16);
py::tuple entry_flag(2); entry_flag[0] = "flag"; entry_flag[1] = list_flag; 
results[0] = entry_flag;
py::tuple entry_ohdu(2); entry_ohdu[0] = "ohdu"; entry_ohdu[1] = list_ohdu; 
results[1] = entry_ohdu;
py::tuple entry_runID(2); entry_runID[0] = "runID"; entry_runID[1] = list_runID; 
results[2] = entry_runID;
py::tuple entry_E0(2); entry_E0[0] = "E0"; entry_E0[1] = list_E0; 
results[3] = entry_E0;
py::tuple entry_n0(2); entry_n0[0] = "n0"; entry_n0[1] = list_n0; 
results[4] = entry_n0;
py::tuple entry_nSavedPix(2); entry_nSavedPix[0] = "nSavedPix"; entry_nSavedPix[1] = list_nSavedPix;
results[5] = entry_nSavedPix;
py::tuple entry_xPix(2); entry_xPix[0] = "xPix"; entry_xPix[1] = list_xPix;
results[6] = entry_xPix;
py::tuple entry_yPix(2); entry_yPix[0] = "yPix"; entry_yPix[1] = list_yPix;
results[7] = entry_yPix;
py::tuple entry_ePix(2); entry_ePix[0] = "ePix"; entry_ePix[1] = list_ePix;
results[8] = entry_ePix;
py::tuple entry_level(2); entry_level[0] = "level"; entry_level[1] = list_level;
results[9] = entry_level;
py::tuple entry_T(2); entry_T[0] = "T"; entry_T[1] = list_T;
results[10] = entry_T;
py::tuple entry_L(2); entry_L[0] = "L"; entry_L[1] = list_L;
results[11] = entry_L;
py::tuple entry_stdT(2); entry_stdT[0] = "stdT"; entry_stdT[1] = list_stdT;
results[12] = entry_stdT;
py::tuple entry_stdL(2); entry_stdL[0] = "stdL"; entry_stdL[1] = list_stdL;
results[13] = entry_stdL;
py::tuple entry_errT(2); entry_errT[0] = "errT"; entry_errT[1] = list_errT;
results[14] = entry_errT;
py::tuple entry_errL(2); entry_errL[0] = "errL"; entry_errL[1] = list_errL;
results[15] = entry_errL;
f.Close();
return_val = results;
            """%selection
        #print code
        partial_result = scipy.weave.inline( code, [ 'fname' ], 
                        headers=['"TFile.h"', '"TTree.h"', '"TObject.h"', '"TObjArray.h"'], 
                        libraries=['Core'],
                        include_dirs=['/usr/include/root/'],
                        library_dirs=['/usr/lib/x86_64-linux-gnu/root5.34/'],
                        #compiler= 'gcc -O3'
                        extra_compile_args=['-O3'],
                        )
        del t
        return dict(partial_result)
    
    def readTree():
        pass
    
    
    def readDict( self, fname ):
        t = timer()
        print 'parseTree'
        if not os.path.exists( fname ):
            print "dict not found", fname
            return None
        if not '.dict' in fname:
            print "file is not dict-file", fname
            return None
        #result = eval(open(fname).read())
        result = np.load(fname).items()
        t()
        #print result['hitSumm']['xPix']
        print "report"
        print len(result['E0'])
        print len(result['nSavedPix'])
        print len(result['xPix'])
        print len(result['ePix'])

        return result
        #self.tree = f.Get("hitSumm")
        #t1 = timer("parseTree")
        #print "parse tree"
        #ROOT.gROOT.ProcessLine("\
            #struct event_t {\
            #Int_t xPix[100000];\
            #Int_t yPix[100000];\
            #Int_t level[100000];\
            #Float_t ePix[100000];\
            #Int_t flag;\
            #Int_t runID;\
            #Int_t ohdu;\
            #Float_t E0;\
            #Float_t n0;\
            #Int_t nSavedPix;\
            #Float_t T;\
            #Float_t L;\
            #Float_t stdT;\
            #Float_t stdL;\
            #Float_t errT;\
            #Float_t errL;\
            #Float_t ecc2;\
            #Float_t ecc3;\
            #Float_t ecc4;\
            #Float_t slope30;\
            #Float_t intp30;\
            #Float_t slopeErr30;\
            #Float_t slope50;\
            #Float_t intp50;\
            #Float_t slopeErr50;\
            #Float_t stdZ10;\
            #Float_t stdZ010;\
            #Float_t stdZerr10;\
            #Float_t stdZ_ex10;\
            #Float_t stdZ0_ex10;\
            #Float_t stdZerr_ex10;\
            #Float_t stdZ30;\
            #Float_t stdZ030;\
            #Float_t stdZerr30;\
            #Float_t stdZ50;\
            #Float_t stdZ050;\
            #Float_t stdZerr50;\
            #Float_t stdZ_ex30;\
            #Float_t stdZ0_ex30;\
            #Float_t stdZerr_ex30;\
            #Float_t stdZ_ex50;\
            #Float_t stdZ0_ex50;\
            #Float_t stdZerr_ex50;\
            #}" )
        #event = ROOT.event_t()
        #self.varsExtra = self.scalarvars + self.listExtra + self.listvars
        #for var in self.varsExtra: 
            #tree.SetBranchAddress( var, AddressOf( event, var ) )
            #self.__dict__[var] = []
            
        #self.tracks = []
        #tracks = self.tracks
        #localvars = self.__dict__
        #t1.factor(text = 'entries', f = 1./tree.GetEntries() )
        
        #for i in xrange( tree.GetEntries() ):
            #trackID = i
            ##if i > 10000: break
            ##t = timer("for")
            #tree.GetEntry(i)
            #if event.nSavedPix > 20: # and event.ohdu == 6:
                #track = Track()
                #track.id = trackID
                #for var in self.scalarvars + self.listExtra:
                    #setattr(track, var, getattr( event, var ) )
                    ##self.__dict__[var] += [ getattr( event, var ) ]
                    ##self.__dict__[var].append( getattr( event, var ) )
                    #localvars[var].append( getattr( event, var ) )
                #for var in self.listvars:
                    #setattr(track, var, [ getattr( event, var )[pix] for pix in xrange( event.nSavedPix ) ] )
                #tracks.append( track )
	#for var in self.varsExtra:
	    #localvars[var] = np.array( localvars[var] )
	#tree.ResetBranchAddress(0)
	##tree.Delete()
	#f.Close()
        #return self


if __name__ == '__main__':
    r = root2numpy()
    r.weaveFull(sys.argv[-1])
    r.weaveFullExtraSelection(sys.argv[-1],'if(nSavedPix<20){continue;}')
    #r.weave2(sys.argv[-1])
    #r.weave(sys.argv[-1])
    r.python(sys.argv[-1])
    #r.readDict( sys.argv[-1] + '.dict.npz' )

