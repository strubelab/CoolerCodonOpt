"""
Make plots for the scores of the optimized sequences after running optimize.py
to generate 25, 50, 100 and 200 sequences
"""

import matplotlib.pyplot as plt
import numpy as np

def get_scores():
    """
    Get the scores from the following files:

    scores_25.txt
    scores_50n.txt
    scores_50nback.txt
    scores.txt
    scores_200n.txt
    scores_200nback.txt

    Example of a file

    Optimized score: -61.5400
    Optimized score: -60.2200
    Optimized score: -61.5200
    Optimized score: -60.6000

    Output
    -----

    Matrix with the scores for each file, with the following columns

    0 - scores_25
    1 - scores_50n
    2 - scores_50nback
    3 - scores
    4 - scores_200n
    5 - scores_200nback

    """

    fnames = [
                '../results/scores/scores_25.txt',
                '../results/scores/scores_50n.txt',
                '../results/scores/scores_50nback.txt',
                '../results/scores/scores.txt',
                '../results/scores/scores_200n.txt',
                '../results/scores/scores_200nback.txt'
                ]

    scores = np.zeros((50,6))

    for i in range(len(fnames)):
        
        with open(fnames[i]) as f:
            lines = f.readlines()
            
            for j in range(len(lines)):
                scores[j,i] = float(lines[j].split()[2])

    return scores


def plot_scores(scores):
    """
    Plot the scores of the tests

    Input
    -----
    scores : numpy array
        Array with the scores for each test
    """

    fig, ax = plt.subplots()

    ax.plot(scores[:,0], label='25 seqs')

    ax.plot(scores[:,1], label='50n seqs')

    ax.plot(scores[:,2], label='50nback seqs')

    ax.plot(scores[:,3], label='100 seqs')

    ax.plot(scores[:,4], label='200n seqs')

    ax.plot(scores[:,5], label='200nback seqs')

    for i in range(6):
        mu = np.mean(scores[:,i])
        sigma = np.std(scores[:,i], ddof=1)

        ax.text(0.7, 0.5+(i*0.05), 
            r'$\mu_%d=%.2f,\ \sigma=%.2f$' % (i, mu, sigma), 
            transform=ax.transAxes, fontsize=15)

        # Middle line
        ax.plot(np.linspace(0,50,10), np.array([mu]*10), linestyle='dotted', 
            color='black')

    ax.legend()

    plt.show()


def plot_50nvsnback(scores):
    """
    Plot the n vs nback versions of the scores

    Input
    -----
    scores : numpy array
        Array with the scores for each test
    """

    fig, ax = plt.subplots()

    ax.plot(scores[:,1], label='50n seqs')

    ax.plot(scores[:,2], label='50nback seqs')

    for i in range(2):
        mu = np.mean(scores[:,i+1])
        sigma = np.std(scores[:,i+1], ddof=1)

        ax.text(0.7, 0.75+(i*0.05), 
            r'$\mu_%d=%.2f,\ \sigma=%.2f$' % (i+1, mu, sigma), 
            transform=ax.transAxes, fontsize=15)

        # Middle line
        ax.plot(np.linspace(0,50,10), np.array([mu]*10), linestyle='dotted', 
            color='black')

    ax.legend()

    plt.show()


def plot_200nvsnback(scores):
    """
    Plot the n vs nback versions of the scores

    Input
    -----
    scores : numpy array
        Array with the scores for each test
    """

    fig, ax = plt.subplots()

    ax.plot(scores[:,4], label='200n seqs')

    ax.plot(scores[:,5], label='200nback seqs')
    
    for i in range(2):
        mu = np.mean(scores[:,i+4])
        sigma = np.std(scores[:,i+4], ddof=1)

        ax.text(0.7, 0.75+(i*0.05), 
            r'$\mu_%d=%.2f,\ \sigma=%.2f$' % (i+4, mu, sigma), 
            transform=ax.transAxes, fontsize=15)

        # Middle line
        ax.plot(np.linspace(0,50,10), np.array([mu]*10), linestyle='dotted', 
            color='black')

    ax.legend()

    plt.show()


def plot_means(scores):
    """
    Plot the means of all the tests
    """

    means = np.array([np.mean(scores[:,i]) for i in range(6)])

    fig, ax = plt.subplots()
    ax.plot(means)

    labels = ['25', '50n', '50nback', '100', '200n', '200nback']

    ax.set_xticks(np.arange(6))
    ax.set_xticklabels(labels)

    for i in range(6):
        ax.text(i, means[i], str(means[i]))

    plt.show()

scores = get_scores()
# plot_scores(scores)
# plot_50nvsnback(scores)
# plot_200nvsnback(scores)
plot_means(scores)
