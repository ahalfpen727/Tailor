if [ -f ~/.envrc ]; then
   source ~/.envrc
fi

if [ -f ~/.condarc ]; then
   source ~/.envrc
fi

set -o notify
export HISTCONTROL=ignoredups

alias h=history
alias less='less -r'
alias whence='type -a'
alias vi='vim'
alias cls='clear'
alias du='du -c -h'
alias du0='du --max-depth=0'
alias du1='du --max-depth=1'
alias df='df -h'
alias more='less'
alias da='date "+%A, %B %d, %Y [%T]"'
alias psg='ps -Af | grep $1'         # requires an argument

#######################
# safety features
#######################
alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
alias ln='ln -i'
alias chown='chown --preserve-root'
alias chmod='chmod --preserve-root'
alias chgrp='chgrp --preserve-root'

#######################
# ls customized
#######################
alias ls='ls -algthrFC --color'
# alias ls='ls -F --color=tty'
alias ld='ls -d */'
alias ll='ls -l'
alias la='ls -A'
alias lla='ls -lA'
alias lr='ls -R'                    # recursive ls
alias lx='ll -BX'                   # sort by extension
alias lz='ll -rS'                   # sort by size
alias lt='ll -rt'                   # sort by date
alias lm='la | more'
alias dir='ls -l --format=long'
alias vdir='ls -l --format=long'
# alias dir='ls --color=auto --format=vertical'
# alias vdir='ls --color=auto --format=long'


#######################
# cd aliases
#######################
alias c:='cd /cygdrive/c'
alias ..="cd .."
alias ..2="cd ../.."
alias ..3="cd ../../.."
alias ..4="cd ../../../.."
alias ..5="cd ../../../../.."

#source ~/.titlebarrc

#######################
# grep customized
#######################
alias grep='grep --color'
alias grep='grep --color=always'
alias egrep='egrep --color=always'
alias fgrep='fgrep --color=always'
#export GREP_OPTIONS='--color=auto' GREP_COLOR='1;32'
export GREP_OPTIONS='--color=always'
export GREP_COLORS="ms=31:mc=31:sl=33:cx=:fn=35:ln=32:bn=32:se=36"
#export GREP_COLORS="ms=01;31:mc=01;31:sl=01;33:cx=:fn=35:ln=32:bn=32:se=36"
#export GREP_COLORS=""

#######################
# functions
#######################

# A function to pipe any command to less:
function so {
eval "$@" |less -I~
}

# extract <file1> <file2> ...
extract() {
    local c e i

    (($#)) || return

    for i; do
        c=''
        e=1

        if [[ ! -r $i ]]; then
            echo "$0: file is unreadable: \`$i'" >&2
            continue
        fi

        case $i in
        #*.t@(gz|lz|xz|b@(2|z?(2))|a@(z|r?(.@(Z|bz?(2)|gz|lzma|xz)))))
        #       c='bsdtar xvf';;
        *.7z)  c='7z x';;
        *.Z)   c='uncompress';;
        *.bz2) c='bunzip2';;
        *.exe) c='cabextract';;
        *.gz)  c='gunzip';;
        *.rar) c='unrar x';;
        *.xz)  c='unxz';;
        *.zip) c='unzip';;
        *)     echo "$0: unrecognized file extension: \`$i'" >&2
               continue;;
        esac

        command $c "$i"
        e=$?
    done

    return $e
}

###############################################################################
# COLOR STUFF
###############################################################################

eval `dircolors -b ~/.dircolors`

# Reset
Color_Off='\e[0m'       # Text Reset

# Regular Colors        "Real Color"
Black='\e[0;30m'        # Black
Red='\e[0;31m'          # Red
Green='\e[0;32m'        # Green
Brown='\e[0;33m'        # Brown
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan
Grey='\e[0;37m'         # Light Grey

# Bold and/or High Intensity
Bold='\e[0m'            # Bold and/or High Intensity
BBlack='\e[1;30m'       # Bold Black OR Dark Grey (Light Black)
BRed='\e[1;31m'         # Bold Red OR Light Red
BGreen='\e[1;32m'       # Bold Green OR Light Green
BBrown='\e[1;33m'       # Bold Brown OR Yellow (Light Brown)
BBlue='\e[1;34m'        # Bold Blue OR Light Blue
BPurple='\e[1;35m'      # Bold Purple OR Light Purple
BCyan='\e[1;36m'        # Bold Cyan OR Light Cyan
BGrey='\e[1;37m'        # Bold Grey OR White (Light Light Grey)

# Underline
UBlack='\e[4;30m'       # Black
URed='\e[4;31m'         # Red
UGreen='\e[4;32m'       # Green
UBrown='\e[4;33m'       # Brown
UBlue='\e[4;34m'        # Blue
UPurple='\e[4;35m'      # Purple
UCyan='\e[4;36m'        # Cyan
UGrey='\e[4;37m'        # Light Grey

# Background
On_Black='\e[40m'       # Black
On_Red='\e[41m'         # Red
On_Green='\e[42m'       # Green
On_Brown='\e[43m'       # Brown
On_Blue='\e[44m'        # Blue
On_Purple='\e[45m'      # Purple
On_Cyan='\e[46m'        # Cyan
On_Grey='\e[47m'        # Light Grey

# High Intensty (Lighter)  (from AIXTERM, emacs does not currently support these)
IBlack='\e[0;90m'       # Dark Grey (Light Black)
IRed='\e[0;91m'         # Light Red
IGreen='\e[0;92m'       # Light Green
IBrown='\e[0;93m'       # Yellow (Light Brown)
IBlue='\e[0;94m'        # Light Blue
IPurple='\e[0;95m'      # Light Purple
ICyan='\e[0;96m'        # Light Cyan
IGrey='\e[0;97m'        # White (Light Light Grey)

# Bold High Intensty  (from AIXTERM, emacs does not currently support these)
BIBlack='\e[1;90m'      # Bold Dark Grey (Light Black)
BIRed='\e[1;91m'        # Bold Light Red
BIGreen='\e[1;92m'      # Bold Light Green
BIBrown='\e[1;93m'      # Bold Yellow (Light Brown)
BIBlue='\e[1;94m'       # Bold Light Blue
BIPurple='\e[1;95m'     # Bold Light Purple
BICyan='\e[1;96m'       # Bold Light Cyan
BIGrey='\e[1;97m'       # Bold White (Light Light Grey)

# High Intensty backgrounds  (from AIXTERM, emacs does not currently support these)
On_IBlack='\e[100m'   # Dark Grey (Light Black)
On_IRed='\e[101m'     # Light Red
On_IGreen='\e[102m'   # Light Green
On_IBrown='\e[103m'   # Yellow (Light Brown)
On_IBlue='\e[104m'    # Light Blue
On_IPurple='\e[105m'  # Light Purple
On_ICyan='\e[106m'    # Light Cyan
On_IGrey='\e[107m'    # White (Light Light Grey)

# set a nice looking prompt:
PS1="\[$BCyan\]\[$On_Purple\]\w\[$Green\]\[$On_Black\] \$\[$BGrey\]"

# Environment  variables
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/pkg/jdk/1.7.0_25/jre/lib/amd64/server
PATH="/home/tr29b/perl5/bin${PATH+:}${PATH}"
export PATH
export PATH=$PATH:/project/umb_triley/cpct/rna-seq/tailor
PERL5LIB="/home/tr29b/perl5/lib/perl5${PERL5LIB+:}${PERL5LIB}"
export PERL5LIB
PERL_LOCAL_LIB_ROOT="/home/tr29b/perl5${PERL_LOCAL_LIB_ROOT+:}${PERL_LOCAL_LIB_ROOT}"
export PERL_LOCAL_LIB_ROOT
PERL_MB_OPT="--install_base \"/home/tr29b/perl5\""
export PERL_MB_OPT
PERL_MM_OPT="INSTALL_BASE=/home/tr29b/perl5"
export PERL_MM_OPT
