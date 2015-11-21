# CS 75 Project: Profile Hidden Markov Models

- Model code located inside `HMMModel.py`
- Viterbi algorithm (standard version) located inside `viterbi.py`
- Viterbi algorithm (diagonalized, wavefront version) located inside `viterbi_diag.py`
- Viterbi algorithm (diagonalized, `with threads`) located inside `viterbi_concur.py`
- Test cases located in the named subdirectory


## Project Abstract
A Viterbi-based algorithm for profile-based alignment is used to accurately determine the optimal alignment between related sequences using position-specific scoring. This project implements a Hidden Markov Model that relies on a Viterbi-based algorithm to test whether specific query sequences belong to a protein family. We apply our implementation within case analyses in which we firstly simulate the process of taking a sequence of unknown origins and searching for a family where said sequence has a domain. Secondly, we use cross-validation to display the relationships between families of the same clan and different clans. Additionally, this project explores the implementation of a parallelized “wave-front” pattern in order to enhance runtime performance through a multi-threaded approach. This pattern is applied by reading each diagonal in our Viterbi-based model as an iteration in the parallelized model. Although parallelizing will not result in changes in asymptotic behavior, parallelizing can make a noticeable difference in actual runtime in practice.


## Acquiring Our Repo Locally
1. Navigate to your terminal and pick a directory
2. `git clone https://github.com/cs75/project.git` inside chosen directory
3. Directory is now available locally on your machine



------

(Below relevant to team members only)

## Instructions on how to clone & sync so you don't need to keep re-downloading from Google Drive


## To get the files
1. Open up your `Terminal` application
2. Type `cd Desktop` and hit enter
3. Copy and paste `git clone https://github.com/cs75/project.git` into your Terminal and hit enter
4. Now the project is on your desktop, and synced up

## To update your local files (when others make changes)
1. The clone above clones the project locally, but it is linked to this remote on GitHub.
2. `git pull origin master`

Alternatively you can just reclone...

## To make changes and push them to remote
1. Change your files
2. In the terminal:
  - `git diff` (make sure changes are correct)
  - `git add .`
  - `git commit -m "i made some changes`
  - `git push origin master`
