import pandas as pd
from utils import getRsids
import glob
from datetime import date

today = date.today()

variantPositions_CTL = pd.read_csv(glob.glob("data/2023-01-10_AIresCTL_genes.csv")[0])
variantRsids_CTL = getRsids(variantPositions_CTL)

variantsAll_CTL = pd.concat([variantPositions_CTL, variantRsids_CTL], axis = 1)
variantsAll_CTL.to_csv('data/' + today.strftime("%Y-%m-%d") + '_AIresCTL_rsids.csv', index = False)

variantPositions_FNF = pd.read_csv(glob.glob("data/2023-01-10_AIresFNF_genes.csv")[0])
variantRsids_FNF = getRsids(variantPositions_FNF)

variantsAll_FNF = pd.concat([variantPositions_FNF, variantRsids_FNF], axis = 1)
variantsAll_FNF.to_csv('data/' + today.strftime("%Y-%m-%d") + '_AIresFNF_rsids.csv', index = False)
