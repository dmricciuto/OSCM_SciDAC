import sys
import argparse
import model_sELM as models

parser = argparse.ArgumentParser()
parser.add_argument('coords', type=float, nargs='*',
                    help="coordinates (lon, lat)")
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
    site_coords = [-72.17, 42.54]
else:
    site_coords = [args.coords[0], args.coords[1]]

if args.pnames_file is not None:
    with open(args.pnames_file) as f:
        pnames_ens = f.read().splitlines()
else:
    pnames_ens = ['nue_tree', 'gdd_crit']  # 'leafcn', 'slatop_decid']

ndim = len(pnames_ens)
nens = args.nens

out_prefix = args.out_prefix

# create model object
model = models.MyModel()

# Desired outputs
#myoutvars = ['gpp', 'lai', 'npp']
print('Generating ensemble')
model.generate_ensemble(nens, pnames_ens, fname=args.psam_file, normalized=False)
assert(model.parm_ensemble.shape[0] == nens)
assert(model.parm_ensemble.shape[1] == ndim)

# Daily output
# model.run_selm(spinup_cycles=4,
#                lon_bounds=[-82, -80], lat_bounds=[40, 42],
#                prefix='regional', pft=0, use_MPI=False)

# Monthly output
# model.run_selm(spinup_cycles=4,
#                lon_bounds=[-85,-80], lat_bounds=[40,45],
#                prefix='regional', do_monthly_output=True,
#                pft=0, use_MPI=True)

# if len(sys.argv) == 1:
#     site_coords = [-72.17, 42.54]  # -72.17 42.54    b'US-Ha1'
# else:
#     assert(len(sys.argv) == 3)
#     site_coords = [float(sys.argv[1]), float(sys.argv[2])]

print("Simulation requested for site coords: ", site_coords)

#model.run_selm(spinup_cycles=0, lon_bounds=[-84.5,-84.2], lat_bounds=[45,45.5], prefix='regional', do_monthly_output=True, deciduous=True)
#model.run_selm(spinup_cycles=4, lon_bounds=[-72.1,-71.9], lat_bounds=[42.3,42.6], prefix='regional', do_monthly_output=True, deciduous=True)
#model.run_selm(spinup_cycles=4, lon_bounds=[-86.6,-86.1], lat_bounds=[39.1,39.6], prefix='regional', do_monthly_output=True) #, deciduous=True)

dellon = 0.1
dellat = 0.1
model.run_selm(spinup_cycles=4,
               lon_bounds=[site_coords[0] - dellon, site_coords[0] + dellon],
               lat_bounds=[site_coords[1] - dellat, site_coords[1] + dellat],
               prefix=out_prefix, do_monthly_output=True,
               pft=0, use_MPI=False)
