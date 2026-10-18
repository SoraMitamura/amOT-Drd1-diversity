"""Rerun frozen GO selections without changing archived results."""
from pathlib import Path
import argparse,subprocess,sys
import pandas as pd
from scipy.stats import false_discovery_control
def main():
 p=argparse.ArgumentParser();p.add_argument('--input-dir',required=True);p.add_argument('--output-dir',default='rerun/GO');p.add_argument('--resamples',type=int,default=9999);p.add_argument('--seed',type=int,default=20261017);a=p.parse_args()
 root=Path(__file__).resolve().parent;out=Path(a.output_dir);tables=[]
 for s in sorted((root/'results/GO').glob('*/selection.json')):
  dest=out/s.parent.name
  subprocess.run([sys.executable,str(root/'go_marker_excluded/analysis.py'),'--input-dir',a.input_dir,'--output-dir',str(dest),'--selection',str(s),'--resamples',str(a.resamples),'--seed',str(a.seed)],check=True)
  t=pd.read_csv(dest/'test_results.csv');t.insert(0,'category',s.parent.name);tables.append(t)
 t=pd.concat(tables,ignore_index=True);t['q_BH_48_tests']=false_discovery_control(t.p_unadjusted.to_numpy());t.to_csv(out/'all_48_tests.csv',index=False)
if __name__=='__main__':main()
