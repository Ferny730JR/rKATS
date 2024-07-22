file(REMOVE_RECURSE
  ".0.9.0"
  "libkatss.0.9.0.dylib"
  "libkatss.dylib"
  "libkatss.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/katss_shared.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
