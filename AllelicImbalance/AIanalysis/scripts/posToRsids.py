import pandas as pd
from utils import getRsids

variantPositions = pd.read_csv("data/AIresCTL.csv")
variantRsids = getRsids(variantPositions)

variantsAll = pd.concat([variantPositions, variantRsids], axis = 1)
variantsAll.to_csv('data/AIresCTL_rsids.csv')
