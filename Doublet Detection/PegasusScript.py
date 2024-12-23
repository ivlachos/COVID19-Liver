

import pegasus as pg
data = pg.read_input(
  '/data/work/Projects/BrScRNAseq/data/BroadTerra/liver.zarr.zip'
)

pg.infer_doublets(data) 

data.obs.to_csv(r'/data/work/Projects/BrScRNAseq/data/BroadTerra/liver.Pegasus.csv')

