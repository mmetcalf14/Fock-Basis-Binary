#!/usr/bin/env bash
rm nohup.out

git init
git add .

echo Which branch are you working in\?
git branch
read Branch

echo Enter your commit message\, please.
read CommitMessage


git commit -m "$CommitMessage"
git push -u origin $Branch
