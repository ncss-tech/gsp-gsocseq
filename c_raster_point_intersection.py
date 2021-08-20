# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 10:23:40 2021

@author: Charles.Ferguson
"""
import arcpy

ws = r'D:\geodata\masks'
pts = r'D:\geodata\fin_project_data4\fin_project_data4\gsp-gsocseq_pe\gsocseq_maps.gpkg\main.gsocseq_maps'
arcpy.management.MakeFeatureLayer(pts, "pts")
arcpy.env.workspace = ws

rasters = arcpy.ListRasters("*")
arcpy.env.overwriteOutput = True
for raster in rasters:
    x = arcpy.Describe(raster).Extent.polygon
    out = raster.replace("mask.tif", "extent.shp")
    arcpy.management.CopyFeatures(x, out)
    arcpy.management.MakeFeatureLayer(out, "poly")
    arcpy.management.SelectLayerByLocation("pts", "INTERSECT", "poly")
    cnt = arcpy.management.GetCount("pts")
    print('{} extent intersects {} points'.format(raster, cnt[0]))
    arcpy.management.Delete("poly")
    

