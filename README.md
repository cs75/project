# project

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
