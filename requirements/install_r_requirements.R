
# --- helper function ------
install.with.check <- function(package) {
  print(sprintf('Package: %s', package))
  if (package == "ADNIMEERGE") {
    print('  - see here for installing ADNIMERGE: https://adni.bitbucket.io/')
  } else if (! package %in% rownames(installed.packages())) {
    print('    - Installing...')
    install.packages(package)
    print('    - Done.')
  } else {
    print('    - Already installed!')
  }
}

# --- set wd -------

install.with.check('this.path')

library(this.path)
setwd(this.dir())

REQUIREMENTS = 'requirements_r.txt'
tbl <- read.table(REQUIREMENTS, header=F)
packages <- tbl$V1

for (package in packages) {
  install.with.check(package)
}
