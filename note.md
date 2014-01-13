PASA installation
=================

MySQL setup
-----------
read-only account = pasa_guest@localhost, password = 'rnaseqasm'

all privilages account = pasa@localhost, password = 'rnaseqasm'

command
-------
create user pasa_guest@localhost identified by 'rnaseqasm';

create user pasa_guest@localhost identified by 'rnaseqasm';

grant select on *.* to pasa_guest@localhost;

grant all on *.* to pasa@localhost;

Perl MySQL module installation
------------------------------

In CPAN,

install DBD::mysql

FASTA3 Suite installation
-------------------------

cd /mnt/data/source

wget http://faculty.virginia.edu/wrpearson/fasta/fasta3/CURRENT.tar.gz

tar xvfz CURRENT.tar.gz

cd fasta-3/src/

make -f ../make/Makefile.linux_sse2 all

cd fasta-3/bin

ln -fs fasta35 fasta

GMAP installation
-----------------

cd /mnt/data/source

wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2013-11-27.tar.gz

tar xvfz gmap-gsnap-2013-11-27.tar.gz

cd gmap-gsnap-2013-11-27/

./configure && make && make install

BLAT installation
-----------------

make sure BLAT in in PATH.
