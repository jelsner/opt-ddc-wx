#!/usr/bin/env python3
"""Download monthly ERA5 hourly weather grids for the conterminous US."""

import argparse
import os


CONUS_AREA = [49.5, -125.0, 24.0, -66.5]  # north, west, south, east
ERA5_VARIABLES = [
    "2m_temperature",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Download monthly ERA5 single-level hourly grids for CONUS. "
            "Requires a configured CDS API key."
        )
    )
    parser.add_argument("--year", type=int, required=True, help="Year to download, e.g. 2020.")
    parser.add_argument(
        "--output-dir",
        default="data/output/conus",
        help="Directory for NetCDF files. Defaults to data/output/conus.",
    )
    parser.add_argument(
        "--months",
        nargs="+",
        type=int,
        default=range(1, 13),
        help="One or more month numbers. Defaults to all months.",
    )
    parser.add_argument(
        "--area",
        nargs=4,
        type=float,
        metavar=("NORTH", "WEST", "SOUTH", "EAST"),
        default=CONUS_AREA,
        help="ERA5 area bounding box. Defaults to CONUS.",
    )
    parser.add_argument(
        "--grid",
        nargs=2,
        type=float,
        metavar=("LAT_RES", "LON_RES"),
        default=None,
        help="Optional output grid spacing in degrees, e.g. --grid 0.5 0.5.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Download even when the destination file already exists.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        import cdsapi
    except ImportError as exc:
        raise SystemExit(
            "Missing Python package 'cdsapi'. Install it and configure your CDS API key "
            "before downloading ERA5 data."
        ) from exc

    client = cdsapi.Client()
    for month in args.months:
        month_str = f"{month:02d}"
        outfile = os.path.join(args.output_dir, f"conus_era5_{args.year}_{month_str}.nc")
        if os.path.exists(outfile) and not args.overwrite:
            print(f"Skipping existing file: {outfile}")
            continue

        request = {
            "variable": ERA5_VARIABLES,
            "product_type": "reanalysis",
            "year": str(args.year),
            "month": [month_str],
            "day": [f"{day:02d}" for day in range(1, 32)],
            "time": [f"{hour:02d}:00" for hour in range(24)],
            "area": args.area,
            "format": "netcdf",
        }
        if args.grid is not None:
            request["grid"] = args.grid

        print(f"Requesting ERA5 CONUS data for {args.year}-{month_str} -> {outfile}")
        client.retrieve("reanalysis-era5-single-levels", request, outfile)


if __name__ == "__main__":
    main()
