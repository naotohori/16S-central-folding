grep '^bond'     16SCD_ninfo.ninfo > bond.ninfo
grep '^angl'     16SCD_ninfo.ninfo > angl.ninfo
grep '^bs-dist'  16SCD_ninfo.ninfo > bs-dist.ninfo
grep '^bs-dihd'  16SCD_ninfo.ninfo > bs-dihd.ninfo
grep '^tbs-dist' 16SCD_ninfo.ninfo > tbs-dist.ninfo
grep '^tbs-angl' 16SCD_ninfo.ninfo > tbs-angl.ninfo
grep '^tbs-dihd' 16SCD_ninfo.ninfo > tbs-dihd.ninfo
grep '^hb-dist'  16SCD_ninfo.ninfo > hb-dist.ninfo
grep '^hb-angl'  16SCD_ninfo.ninfo > hb-angl.ninfo
grep '^hb-dihd'  16SCD_ninfo.ninfo > hb-dihd.ninfo

# Edit tbs-dist.ninfo to add +/- from 16SCD_16SCD_unprocessed_stacks.dat

# Edit hb-dist.ninfo to add nHB and atom names from output of hb_list_convert.py (16SCD.hb.list)

cat  bond.ninfo \
     angl.ninfo \
     bs-dist.ninfo \
     bs-dihd.ninfo \
     tbs-dist.edit.ninfo \
     tbs-angl.ninfo \
     tbs-dihd.ninfo \
     hb-dist.edit.ninfo \
     hb-angl.ninfo \
     hb-dihd.ninfo \
     > 16SCD.ninfo
