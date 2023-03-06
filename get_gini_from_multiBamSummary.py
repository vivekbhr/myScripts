#!/usr/bin/env python
import sys
sys.path.append("/hpc/hub_oudenaarden/vbhardwaj/programs/sincei/sincei")
from Utilities import gini
import pandas as pd
from scipy import sparse, io
import re

tab = pd.read_csv(sys.argv[1], sep="\t")
tab.columns = [re.sub("'", '', x) for x in tab.columns]
varnames = tab['#chr'].map(str) + '_' + tab['start'].map(str) + '_' + tab['end'].map(str)
tab2 = tab.iloc[:, 3:(tab.shape[1])].transpose()
tab2.columns = varnames
tab3=sparse.csr_matrix(tab2)
gini_list = [gini(i, tab3) for i in range(tab3.shape[0])]
out = pd.DataFrame({'sample': tab2.index.to_list(),
                   'gini': gini_list})
out.to_csv(sys.argv[2], sep="\t", index=False)
