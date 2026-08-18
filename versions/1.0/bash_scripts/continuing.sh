#!/bin/bash
read -p "Check quality before continuing. \n did you check? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "continue"
fi
