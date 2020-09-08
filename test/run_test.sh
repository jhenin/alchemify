
echo "------------ Running NAMD PSF test ------------------"
../alchemify  poea.psf out.psf poea.pdb | tee out.log
diff out.psf ref_out.psf > out.diff

if [ $? -ne 0 ]
then
  echo "Output of NAMD PSF test differs from reference. Please check out.diff and out.log".
fi

echo
echo "------------ Running CHARMM EXT PSF test ------------------"
../alchemify  ext.psf  out_ext.psf ext.pdb | tee out_ext.log
diff out_ext.psf ref_out_ext.psf > out_ext.diff

if [ $? -ne 0 ]
then
  echo "Output of CHARMM EXT PSF test differs from reference. Please check out_ext.diff and out_ext.log".
fi

echo
echo "Done. Run rm -f out* to clean up output files."
