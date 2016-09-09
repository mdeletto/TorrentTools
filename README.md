# TorrentTools

TorrentTools is a suite of bioinformatics related tools for processing IonTorrent data.  At this time, the suite includes:

1. Uniformity - a TorrentSuite plugin that calculate uniformity per pool.  Useful for determining library issues specific to one or more primer pools.
2. IR_download.py - a tool to download .zip files from an IonReporter server, given an analysis name and optional ID.
3. igv_snapshot_controller.py - a tool to take batch snapshots from IGV.  Requires a local IGV instance already be running.
4. tumor_normal_coverage_analysis.py - a tool for calculating coverage statistics, as well as target region dropouts in matched tumor and normal samples.
