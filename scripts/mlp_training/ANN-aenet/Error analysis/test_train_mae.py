"""
This script reads the 'train.out' file, to extract and plot the Mean Absolute Error (MAE) values 
for both the training and test sets across different iterations of ML potential training. 
"""

import numpy
import matplotlib.pyplot as plt

def get_train_test_MAE():
    '''
        get_train_test_MAE
        Read train.out, extract MAE for training and test set.
    '''
    ffile = open('train.out','r')
    lines = ffile.readlines()
    ffile.close()
    #
    mae_train = []
    mae_test  = []
    step = []
    # Total energy
    for a,line in enumerate(lines):
        splitt = line.split()
        if len(splitt) > 1:
            if splitt[0] == '|------------TRAIN-----------|' and splitt[1] == '|------------TEST------------|': 
                i = 252 # choose to skip the initial steps where the error drops sharply
                stop = False
                while not stop:
                    i += 1
                    splitt2 = lines[a+i].split()
                    if len(splitt2) == 0:
                        stop = True
                    else:
                        mae_train.append(float(splitt2[1]))
                        mae_test.append(float(splitt2[3]))
                        step.append(i)
    return numpy.array(mae_train), numpy.array(mae_test), numpy.array(step)

def plot_stuff():
    # Plot training and test MAE
    mae_train, mae_test, step = get_train_test_MAE()
    plt.plot(step, mae_train,'-',color='green',linewidth=2.0,label="Training set")
    plt.plot(step, mae_test,'-',color='red',linewidth=2.0,label="Test set")
    plt.xlabel('Training iteration',fontsize=15,labelpad=10)
    plt.ylabel('MAE [eV]',fontsize=15,labelpad=10)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.legend(fontsize=15)
    plt.tight_layout()
    plt.savefig('test_train_MAE.png')

def main():
    plot_stuff()


if __name__ == '__main__':
    main()

