# iTerm
(superior terminal)
https://iterm2.com

## Change from zsh to bash
chsh -s /bin/bash

# ssh
Copy over existing keys
```
rsync -Prav user@oldComputer:~/.ssh ~/Desktop
mv ~/.ssh ~/.sshORIGINAL # keep existing just in case
mv ~/Desktop/.ssh ~/.ssh
ssh-add -K # enter your ssh passphrase
```

# brew
OS X package manager for installing command line (and sometimes GUI) applications
https://brew.sh
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

## Create ~/.bash_profile and add the following
```
export BASH_SILENCE_DEPRECATION_WARNING=1 # turns off zsh nag
export PATH="/opt/homebrew/bin/:$PATH" # puts brew into the used path

# append history after each invocation to maintain history across windows
export PROMPT_COMMAND="history -a; history -c; history -r; $PROMPT_COMMAND"
export HISTTIMEFORMAT="%d/%m/%y %T "
```

## KeepingYouAwake
(caffeinate menu bar toggle)
```
brew install --cask keepingyouawake
```

## Miniconda
```
brew install --cask miniconda
conda init bash
```

# R
Install pkg (make sure it is Apple silicon arm64) from https://cran.r-project.org/bin/macosx/

## rstudio desktop
https://posit.co/download/rstudio-desktop/

# energiza
https://appgineers.de/energiza/
Allows you to keep battery state of charge within user-set limits (keeping battery regularly at 100% or below 20% is damaging for the battery)
