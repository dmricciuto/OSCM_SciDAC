#!/usr/bin/env python

import argparse

import utils
import model_sELM as selm


parser = argparse.ArgumentParser()
parser.add_argument('coords', type=float, nargs='*',
                    help="coordinates (lon, lat)")
parser.add_argument("-t", "--site_name", dest="site_name", type=str, default='US-UMB',
                    help="Site name")
parser.add_argument("-p", "--pnames_file", dest="pnames_file", type=str, default=None,
                    help="Parameter names file")
parser.add_argument("-s", "--psam_file", dest="psam_file", type=str, default='',
                    help="Parameter samples file")
parser.add_argument("-n", "--nens", type=int, dest="nens", default=3,
                    help="Ensemble size")
parser.add_argument("-o", "--out_prefix", dest="out_prefix", type=str, default='regional',
                    help="Output file prefix")
args = parser.parse_args()

if len(args.coords) == 0:
    lon, lat = utils.get_lonlat(args.site_name)
    site_coords = [lon, lat]
else:
    site_coords = [args.coords[0], args.coords[1]]
print("Simulation requested for site coords: ", site_coords)

if args.pnames_file is not None:
    with open(args.pnames_file) as f:
        pnames_ens = f.read().splitlines()
else:
    pnames_ens = ['nue_tree', 'gdd_crit']  # 'leafcn', 'slatop_decid']
print("List of perturbed parameters: ", pnames_ens)

ndim = len(pnames_ens)
nens = args.nens

out_prefix = args.out_prefix

# create model object
model = selm.MyModel()


print('Generating ensemble')
model.generate_ensemble(nens, pnames_ens, fname=args.psam_file, normalized=False)
assert(model.parm_ensemble.shape[0] == nens)
assert(model.parm_ensemble.shape[1] == ndim)





dellon = 0.1
dellat = 0.1
model.run_selm(spinup_cycles=4,
               lon_bounds=[site_coords[0] - dellon, site_coords[0] + dellon],
               lat_bounds=[site_coords[1] - dellat, site_coords[1] + dellat],
               prefix=out_prefix, do_monthly_output=False,
               pft=-1, use_MPI=False)
