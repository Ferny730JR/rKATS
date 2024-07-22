file(REMOVE_RECURSE
  "libkatss.a"
  "libkatss.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/katss_static.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
