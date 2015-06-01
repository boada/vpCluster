rm -rf redux/*.log redux/VP_*
rm -rf data/flat data/trace
rm -f data/*/*.{fits,fmod,pmod,dist} data/*/detect.*
find . -name "*.log" -exec rm {} \;
