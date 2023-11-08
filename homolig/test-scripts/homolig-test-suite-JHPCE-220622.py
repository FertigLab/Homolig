#Shell wrapper to run test of Homolig algorithm 
# 14 June 2022 AAG 
#$ -cwd
#$ -m e
#$ -M agirgis3@jhmi.edu
#$ -l mem_free=4G,h_vmem=4G,h_fsize=150G
#$ -pe local 2
#$ -o /users/agirgis/job-logs/
#$ -e /users/agirgis/job-logs/

import os
import homolig as hg

inputDir = '/users/agirgis/homolig/homolig/data/test_data'

print(os.getcwd())
fail = 0


#Paired TCR test
inputFile = inputDir +  '/tcr_test.csv'
print('Running paired chain TCR test... ----------------')
def main():
      adata = hg.homolig(input_file = inputFile,
        seq_type = 'tcr',
        chains = 'paired',
        metric='aadist',
        species='human')

try:
    if __name__ == '__main__':
        main()
    print('Successfully passed test. -----------------')
except:
    print('---------------- FAILED TEST!!!! ----------------')
    fail += 1
    pass


#Alpha chain tcr test
print('Running alpha chain TCR test... ----------------')
inputFile = inputDir +  '/homolig_filt_2022-06-14_TRA_HCC20.csv'
def main():
      adata = hg.homolig(input_file = inputFile,
        seq_type ='tcr',
        chains ='alpha',
        metric='aadist',
        species='human')

try:
    if __name__ == '__main__':
        main()
    print('Successfully passed test. -----------------')
except:
    print('---------------- FAILED TEST!!!! ----------------')
    fail += 1
    pass


#Heavy chain bcr test
inputFile = inputDir + '/homolig_filt_2022-06-14_IGH_HCC04.csv'
print('Running heavy chain BCR test... ----------------')
def main():
       adata = hg.homolig(input_file = inputFile,
        seq_type = 'bcr',
        chains = 'heavy',
        metric='aadist',
        species='human')

try:
    if __name__ == '__main__':
        main()
    print('Successfully passed test. -----------------')
except:
    print('---------------- FAILED TEST!!!! ----------------')
    fail += 1
    pass

#Light chain bcr test
inputFile = inputDir +  '/homolig_filt_2022-06-14_IGL_HCC04.csv'
print('Running light chain BCR test... -----------------')
def main():
       adata = hg.homolig(input_file = inputFile,
        seq_type = 'bcr',
        chains = 'light',
        metric='aadist',
        species='human')

try:
    if __name__ == '__main__':
        main()
    print('Successfully passed test. -----------------')
except:
    print('---------------- FAILED TEST!!!! ----------------')
    fail += 1
    pass


#Beta chain tcr test
print('Running beta chain TCR test... -----------------')
inputFile = inputDir + '/homolig_filt_2022-06-14_TRB_HCC04.csv'
def main():
       adata = hg.homolig(input_file = inputFile,
        seq_type = 'tcr',
        chains = 'beta',
        metric='aadist',
        species='human')

try:
    if __name__ == '__main__':
        main()
    print('Successfully passed test. -----------------')
except:
    print('---------------- FAILED TEST!!!! ----------------')
    fail += 1
    pass

#Paired BCR test
print('Running paired chain BCR test... ----------------')
inputFile = inputDir +  '/paired_bcr_test.csv'
def main():
       adata = hg.homolig(input_file = inputFile,
        seq_type = 'bcr',
        chains = 'paired',
        metric='aadist',
        species='human')

try:
    if __name__ == '__main__':
        main()
    print('Successfully passed test. -----------------')
except:
    print('---------------- FAILED TEST!!!! ----------------')
    fail += 1
    pass


#Summarize test results
print('Failed ' + str(fail) + ' of 6 total tests.')
