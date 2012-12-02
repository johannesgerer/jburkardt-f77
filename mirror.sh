#!/bin/bash
declare -A n

n[c]=C
n[f77]="FORTRAN 77"
n[java]=Java
n[m]=MATLAB
n[r]=R
n[cpp]=C++
n[f]="FORTRAN 90"
n[math]=Mathematica
n[py]=Python

for i in "${!n[@]}"
do
		s="${i}_src"
		src="http://people.sc.fsu.edu/~jburkardt/$s/$s.html"
		repo="jburkardt-$i"
		
		#curl -u 'johannesgerer' https://api.github.com/user/repos -d @-
		cat <<EOF
{"name":"$repo"
,"description":"An official Git Mirror of John Burkardt's great collection of ${n[$i]} Software"
,"homepage":"$src"
}
EOF

if !(cd "$repo"); then
		git clone git@ghj:johannesgerer/$repo
fi
(cd $repo && (\
#				wget -m -np --cut-dirs=2 -nH $src &&						\
				cp "${s}.html" README.md &&												\
				git add -A &&																		\
				git commit -m "automatic update from $src" &&		\
		    git gc																				 && \
				git push origin master
				)
)
done


