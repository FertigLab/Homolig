#Testing various homolig modes using test data.
import os
import homolig as hg


os.chdir('/home/aag7319/homolig/homolig/data/test_data')
print(os.getcwd())
fail = 0

#Heavy chain bcr test
print('Running heavy chain BCR test... ----------------')
def main():
    adata = hg.homolig(input_file = './homolig_filt_2022-06-14_IGH_HCC04.csv',
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

#Paired TCR test
print('Running paired chain TCR test... ----------------')
def main():
    adata = hg.homolig(input_file = './tcr_test.csv',
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
def main():
    adata = hg.homolig(input_file = './homolig_filt_2022-06-14_TRA_HCC20.csv',
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
print('Running heavy chain BCR test... ----------------')
def main():
    adata = hg.homolig(input_file = './homolig_filt_2022-06-14_IGH_HCC04.csv',
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
print('Running light chain BCR test... -----------------')
def main():
    adata = hg.homolig(input_file = './homolig_filt_2022-06-14_IGL_HCC04.csv',
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
def main():
    adata = hg.homolig(input_file = './homolig_filt_2022-06-14_TRB_HCC04.csv',
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
def main():
    adata = hg.homolig(input_file = './paired_bcr_test.csv',
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