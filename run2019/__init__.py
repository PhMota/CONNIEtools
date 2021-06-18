from iarray import Selection

M_Si_amu = 28.085 #atomic mass unit

runIDexcluded = [6415, 6475, 6499, 6926, 6927, 6031, 6093, 6096, 6139, 7074, 7222, 7226, 7374]
excluded_sel = Selection.prod( runIDexcluded, key="runID", op="!=" )

excluded_hdu = Selection.sum( [2,3,4,5,8,9,13,14], key="hdu", op="==" )
excluded_ohdu = Selection.sum( [2,3,4,5,8,9,13,14], key="ohdu", op="==" )
off = Selection( '(runID>6226) & (runID<6975)' )
on = Selection( '(runID>6030) & (runID<6227)' ) + '(runID>6974) & (runID<7522)'

off_excl = off * excluded_sel
on_excl = on * excluded_sel

data_files = '/share/storage2/connie/DAna/nuCatalogs/shape_*_data_[6-7]*_to_*_v4.0.root'
sim_files = '/share/storage2/connie/DAna/nuCatalogs/draw_all*.root'
match_files = '/share/storage2/connie/DAna/nuCatalogs/match_*_sim_[6-7]*_to_*_v4.0.root'