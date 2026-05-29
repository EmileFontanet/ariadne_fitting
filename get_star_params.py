#!/usr/bin/env python3

import argparse
from pathlib import Path
import os
from astroquery.simbad import Simbad
from astroARIADNE.fitter import Fitter
from astroARIADNE.star import Star
import traceback
import pandas as pd
DEFAULT_MODELS = [
    "phoenix",
    "btsettl",
    "btnextgen",
    "btcond",
    "kurucz",
    "ck04",
]


def query_star(star_name: str):
    simbad = Simbad()
    simbad.add_votable_fields("ids")

    result = simbad.query_object(star_name)

    if result is None:
        raise ValueError(f"Star '{star_name}' not found in SIMBAD")

    ra = result["ra"][0]
    dec = result["dec"][0]

    gaia_id = None
    for identifier in result["ids"][0].split("|"):
        identifier = identifier.strip()
        if identifier.lower().startswith("gaia dr3"):
            gaia_id = int(identifier[8:].strip())
            break

    return ra, dec, gaia_id


def create_fitter(star_name: str, ra: str, dec: str, gaia_id: int | None,
                  output_dir: Path, n_samples: int, prior_table=None):
    star = Star(star_name, ra, dec, g_id=gaia_id)

    engine = 'dynesty'
    nlive = 500
    dlogz = 0.5
    bound = 'multi'
    sample = 'rwalk'
    threads = 8
    dynamic = False

    fitter = Fitter()
    fitter.star = star
    setup = [engine, nlive, dlogz, bound, sample, threads, dynamic]
    fitter.setup = setup
    if prior_table is not None:
        fitter.prior_setup = {
            'teff': ('normal', prior_table.Teff.values[0], prior_table.Teff_err.values[0]),
            'logg': ('normal', prior_table.logg.values[0], prior_table.logg_e.values[0]),
            'z': ('normal', prior_table.feh.values[0], prior_table.feh_e.values[0]),
            'dist': ('default'),
            'rad': ('default'),
            'Av': ('fixed', 0.0)
        }
        print('Set priors from table')
        print(fitter.prior_setup)
    else:
        fitter.prior_setup = {
            'teff': ('default'),
            'logg': ('default'),
            'z': ('default'),
            'dist': ('default'),
            'rad': ('default'),
            # 'Av': ('default')
            'Av': ('fixed', 0.0)
        }
    # sergio_params = pd.read_csv('../results_dr3.rdb', sep = '\t', skiprows=[1])
    # fitter.prior_setup = {
    #         'teff': ('normal', sergio_params.teff.iloc[i], sergio_params['erteff2'].iloc[i]),
    #     'logg': ('normal', sergio_params['logg_dr3'].iloc[i], sergio_params['erlogg_dr3'].iloc[i]),
    #     'z': ('normal', sergio_params['feh'].iloc[i], sergio_params['erfeh2'].iloc[i]),
    #     'dist': ('default'),
    #     'rad': ('default'),
    #     'Av': ('fixed', 0.0)
    #     }
    fitter.av_law = "fitzpatrick"
    fitter.out_folder = str(output_dir)
    fitter.bma = True
    fitter.models = DEFAULT_MODELS
    fitter.n_samples = n_samples

    return fitter


def run_fit(star_name: str, n_samples: int, prior_table=None, output_dir=Path("results")):

    print(f"Querying SIMBAD for {star_name}")
    ra, dec, gaia_id = query_star(star_name)

    fitter = create_fitter(
        star_name=star_name,
        ra=ra,
        dec=dec,
        gaia_id=gaia_id,
        output_dir=output_dir,
        n_samples=n_samples,
        prior_table=prior_table
    )

    fitter.initialize()
    fitter.fit_bma()


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "star",
        help="Star name resolvable by SIMBAD"
    )

    parser.add_argument(
        "-n",
        "--samples",
        type=int,
        default=10000,
        help="Number of posterior samples"
    )

    parser.add_argument(
        "--env",
        type=str,
        default="local",
        help="Environment (local, cluster)"
    )

    return parser.parse_args()


def main():
    args = parse_args()
    params = pd.read_csv('714_params.csv')
    params = params[params.HD == args.star]
    if (args.env == "cluster"):
        output_dir = Path(f"/srv/scratch/fontanee/ariadne_results/{args.star}")
    else:
        output_dir = Path('results/'+args.star)
    best_fit_file = output_dir / "best_fit_average.dat"

    if best_fit_file.exists():
        print(f"{args.star}: fit already exists")
        return

    output_dir.mkdir(exist_ok=True)
    try:
        run_fit(args.star, args.samples,
                prior_table=params, output_dir=output_dir)
    except Exception as exc:
        print(f"Error processing {args.star}: {exc}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
