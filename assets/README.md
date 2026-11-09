# Reference geometry

The fsaverage5_lh and fsaverage5_rh files contain 10,242 vertices and 20,480 faces per hemisphere, with binary cortical masks, for vertex-level maps. These are reference geometry only.

The surface files contain the project's inflated cortical reference meshes and atlas labels: 163,842 vertices and 327,680 triangular faces per hemisphere. Fields are vertices, faces (MATLAB one-based), label, label_ids and label_names. They contain no participant thickness, covariates, fitted maps or statistical results. Atlas files were converted from the surface_500_L_inflated and surface_500_R_inflated reference templates. Match ROI values by hemisphere-prefixed atlas name rather than relying on a fixed 307/308-column split.

Geographic files contain boundary coordinates and region names only. China outlines retain the bou-series national/provincial geometry and separate maritime-line segments used in the project. World and USA files are the project ECharts map geometries converted into plain GeoJSON. No participant locations, city counts or culture scores are included. Coordinates are longitude/latitude; projection is applied by figures/map_geometry.py. Source-file conversion does not confer a new license on third-party geography or reference templates.
