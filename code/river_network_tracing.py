import arcpy

arcpy.env.overwriteOutput = True

nd = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/hybas_lev2_networks/hybas_lev2_networks_ND_as'
source_pts = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/vmap_sources_as'
mouth_pts = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/vmap_mouths_as'
out_lines = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/as_final'
# ord_sources = 'C:/Users/kcher/Desktop/disser/vmap_prepar/vmap_prepar.gdb/ordered_sources'
batch_size = 100
basins = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/hybas_lev2_as'
endo_basins = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/hybas_lev6_nosa_endo'
sinks = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/sinks_snap_sel'
nd_junctions = 'C:/Users/kcher/Desktop/disser/vmap_final/vmap_final.gdb/hybas_lev2_networks/hybas_lev2_networks_ND_as_Junctions'


src_lyr = 'src_lyr'
arcpy.MakeFeatureLayer_management(source_pts, src_lyr)
mth_lyr = 'mth_lyr'
arcpy.MakeFeatureLayer_management(mouth_pts, mth_lyr)
sinks_lyr = 'sinks_lyr'
arcpy.MakeFeatureLayer_management(sinks, sinks_lyr)

endo_buf = 'in_memory/endo_buf'
endo_buf = arcpy.Buffer_analysis(endo_basins, endo_buf, '100 Meters', dissolve_option='ALL')

endo_diss = 'in_memory/endo_diss'
arcpy.Dissolve_management(endo_buf, endo_diss, multi_part='SINGLE_PART')

endo_diss_bound = 'in_memory/endo_diss_bound'
arcpy.PolygonToLine_management(endo_diss, endo_diss_bound)

junc = 'in_memory/junc'
arcpy.Select_analysis(nd_junctions, junc)
junc_lyr = 'junc_lyr'
arcpy.MakeFeatureLayer_management(junc, junc_lyr)

cursor =  arcpy.da.SearchCursor(basins, ['SHAPE@', 'OBJECTID'])
nbasins = len(list(cursor))
arcpy.AddMessage(str(nbasins) + ' basins')
cursor =  arcpy.da.SearchCursor(basins, ['SHAPE@', 'OBJECTID'])

create_out_dataset = True
batch_cnt = 0
for row in cursor:

    bas_id = row[1]
    arcpy.AddMessage("Processing basin " + str(bas_id) + " out of " + str(nbasins))

    arcpy.AddMessage("\tCreating Closest facility layer...")
    nd_lyr = 'nd_lyr'
    arcpy.MakeClosestFacilityLayer_na(nd, nd_lyr, 'Length', travel_from_to='TRAVEL_FROM',
                                      default_number_facilities_to_find=1, hierarchy='USE_HIERARCHY')
    
    basin = 'in_memory/basin'
    arcpy.Select_analysis(basins, basin, 'OBJECTID = ' + str(bas_id))
    
    arcpy.SelectLayerByLocation_management(src_lyr, 'INTERSECT', basin, "", "NEW_SELECTION")
    arcpy.SelectLayerByLocation_management(sinks_lyr, 'INTERSECT', basin, "", "NEW_SELECTION")
    arcpy.SelectLayerByLocation_management(mth_lyr, 'INTERSECT', basin, "5000 Meters", "NEW_SELECTION")
    
    src_sel = 'in_memory/src_sel'
    mth_sel = 'in_memory/mth_sel'
    sinks_sel = 'in_memory/sinks_sel'
    arcpy.Select_analysis(src_lyr, src_sel)
    arcpy.Select_analysis(mth_lyr, mth_sel)
    arcpy.Select_analysis(sinks_lyr, sinks_sel)
    
    arcpy.AddLocations_na(nd_lyr, 'Facilities', mth_sel, "", "", append="CLEAR")
    arcpy.AddLocations_na(nd_lyr, 'Facilities', sinks_sel, "", "")
    arcpy.AddLocations_na(nd_lyr, 'Line Barriers', endo_diss_bound, "", "")

    arcpy.AddMessage("\tOrdering sources...")

    oids_raw = arcpy.da.TableToNumPyArray(src_sel, 'OBJECTID')
    oids = []
    for oid in oids_raw:
        oids.append(oid[0])
    
    n_batches = len(oids) // batch_size
    if len(oids) % batch_size != 0:
        n_batches += 1

    arcpy.AddMessage("\t\tNumber of sources: " + str(len(oids)))
    arcpy.AddMessage("\t\tBatch size: " + str(batch_size))
    arcpy.AddMessage("\t\tNumber of batches: " + str(n_batches))

    create_basin_dataset = True
    valid_sources = 'in_memory/valid_sources'
    for i in range(n_batches):
        arcpy.AddMessage('\t\tProcessing batch ' + str(i+1) + ' out of ' + str(n_batches))
        
        batch_oids = oids[i*batch_size:(i+1)*batch_size]
        
        batch = 'in_memory/batch'
        oids_str = str(tuple(batch_oids))
        if len(batch_oids) < 2:
            oids_str = oids_str.replace(',', '')
        arcpy.Select_analysis(src_sel, batch, 'OBJECTID IN ' + oids_str)

        # arcpy.AddMessage('\t\t Solving...')
        arcpy.AddLocations_na(nd_lyr, 'Incidents', batch, "", "", append="CLEAR")
        try:
            arcpy.Solve_na(nd_lyr)
        except:
            continue

        curr_incidents = 'in_memory/curr_incidents'
        curr_routes = 'in_memory/curr_routes'
        arcpy.Select_analysis('Incidents', curr_incidents)
        arcpy.Select_analysis('Routes', curr_routes)

        arcpy.AddField_management(curr_incidents, 'IncidentID', 'LONG')
        arcpy.CalculateField_management(curr_incidents, 'IncidentID', '[ObjectID]+' + str(batch_size) + '*' + str(i))

        # incident_oids = arcpy.da.TableToNumPyArray(curr_incidents, 'IncidentID')
        # route_iids = arcpy.da.TableToNumPyArray(curr_routes, 'IncidentID')
        # print '\t\t\tIncident IDs: ', list(incident_oids)
        # print '\t\t\tRoute incident IDs: ', list(route_iids)

        arcpy.JoinField_management(curr_incidents, 'IncidentID', curr_routes, 'IncidentID', ['Total_Length'])

        incidents_sel = 'in_memory/incidents_sel'
        arcpy.Select_analysis(curr_incidents, incidents_sel, '"Total_Length" IS NOT NULL')

        oids_sel_raw = arcpy.da.TableToNumPyArray(incidents_sel, 'OBJECTID')
        oids_sel = []
        for oid in oids_sel_raw:
            oids_sel.append(oid[0])

        arcpy.AddMessage('\t\t\tNumber of valid sources: ' + str(len(oids_sel)))

        if create_basin_dataset:
            arcpy.Select_analysis(incidents_sel, valid_sources)
            create_basin_dataset = False
        else:
            arcpy.Append_management(incidents_sel, valid_sources, "NO_TEST")

        arcpy.Delete_management(curr_incidents)
        arcpy.Delete_management(curr_routes)
        arcpy.Delete_management(incidents_sel)

    ordered_sources = 'in_memory/ordered_sources'
    arcpy.Sort_management(valid_sources, ordered_sources, [["Total_Length", "ASCENDING"]])
    # arcpy.Select_analysis(ordered_sources, ord_sources)

    arcpy.AddMessage("\tTracing river network...")

    oids_ord_raw = arcpy.da.TableToNumPyArray(ordered_sources, 'OBJECTID')
    oids_ord = []
    for oid in oids_ord_raw:
        oids_ord.append(oid[0])

    n_batches = len(oids_ord) // batch_size
    if len(oids_ord) % batch_size != 0:
        n_batches += 1

    arcpy.AddMessage("\t\tNumber of valid sources: " + str(len(oids_ord)))
    arcpy.AddMessage("\t\tBatch size: " + str(batch_size))
    arcpy.AddMessage("\t\tNumber of batches: " + str(n_batches))

    routes = 'in_memory/routes'
    create_basin_dataset = True
    for i in range(n_batches):
        arcpy.AddMessage('\t\tProcessing batch ' + str(i + 1) + ' out of ' + str(n_batches))

        batch_oids = oids_ord[i * batch_size:(i + 1) * batch_size]

        batch = 'in_memory/batch'
        oids_str = str(tuple(batch_oids))
        if len(batch_oids) < 2:
            oids_str = oids_str.replace(',', '')
        arcpy.Select_analysis(ordered_sources, batch, 'OBJECTID IN ' + oids_str)

        # arcpy.AddMessage('\t\t Solving...')
        arcpy.AddLocations_na(nd_lyr, 'Incidents', batch, "", "", append="CLEAR")
        try:
            arcpy.Solve_na(nd_lyr)
        except:
            continue
        # arcpy.AddMessage('\t\t Splitting lines...')

        routes_split = 'in_memory/routes_split'
        arcpy.SelectLayerByLocation_management(junc_lyr, 'INTERSECT', 'Routes', selection_type='NEW_SELECTION')
        arcpy.SplitLineAtPoint_management('Routes', junc_lyr, routes_split, '1 Meter')
        arcpy.AddLocations_na(nd_lyr, 'Facilities', junc_lyr, "", "")

        if create_basin_dataset:
            arcpy.Select_analysis(routes_split, routes)
            create_basin_dataset = False
        else:
            arcpy.Append_management(routes_split, routes, "NO_TEST")

    arcpy.AddMessage('\tDeleting Identical features...')
    arcpy.DeleteIdentical_management(routes, 'Shape', "10 Meters")

    arcpy.AddMessage('\tUnsplitting lines...')
    routes_unsplit = 'in_memory/routes_unsplit'
    arcpy.UnsplitLine_management(routes, routes_unsplit)

    arcpy.AddMessage('\tFlipping lines...')
    arcpy.FlipLine_edit(routes_unsplit)

    arcpy.AddMessage('\tFinal planarizing...')
    routes_ftl = 'in_memory/routes_ftl'
    arcpy.FeatureToLine_management(routes_unsplit, routes_ftl)

    arcpy.DeleteIdentical_management(routes_ftl, 'Shape', "10 Meters")

    routes_fix = 'in_memory/routes_fix'
    arcpy.AddMessage('\tFinal unsplitting...')
    arcpy.UnsplitLine_management(routes_ftl, routes_fix)

    arcpy.AddMessage('\tAppending basin routes...')
    if create_out_dataset:
        arcpy.Select_analysis(routes_fix, out_lines)
        create_out_dataset = False
    else:
        arcpy.Append_management(routes_fix, out_lines, "NO_TEST")
    del routes
    del routes_fix

# arcpy.AddMessage('Planarizing lines...')
# arcpy.FeatureToLine_management(out_lines, out_lines_planar, "10 Meters")
#
# arcpy.AddMessage('Deleting Identical features...')
# arcpy.DeleteIdentical_management(out_lines, 'Shape', "10 Meters")

# arcpy.AddMessage('Flipping lines...')
# arcpy.FlipLine_edit(out_lines)
