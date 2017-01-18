scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161119/raw/rawgrnam13.csv .

scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161119/raw/rawm13.csv .

#############

scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161209_miseqforhiseq/rawhi1/raw.csv ./raw20161209_1.csv

scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161209_miseqforhiseq/rawhi2/raw.csv ./raw20161209_2.csv

scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161209_miseqforhiseq/rawhi3/raw.csv ./raw20161209_3.csv


scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161209_miseqforhiseq/all.bed ./grna.bed

#worked, but different way. why did I do it this way?
scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161209_miseqforhiseq/rawhi3/tot.grna .

#new, untested...
scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_20161209_miseqforhiseq/rawhi3/grna.sam tot.grna



scp -o proxycommand="ssh gate.ebi.ac.uk proxy %h" \
mahogny@ebi-001.ebi.ac.uk:/homes/mahogny/common/data/henriksson_jhuma20161222/raw/raw.csv ./henriksson_jhuma20161222.csv

