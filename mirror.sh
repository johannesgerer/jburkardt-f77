#!/bin/bash
declare -A n

n[c]=C
n[m]=MATLAB
n[cpp]=C++
n[f]="FORTRAN 90"
n[py]=Python

n[f77]="FORTRAN 77"
n[r]=R
n[math]=Mathematica
n[java]=Java

echo  "${!n[@]}"
for i in "${!n[@]}" # "$@"
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
				# wget "$src" -O /dev/null \
				# 		&& echo rm \* -rf \
				# 		wget -m -np --cut-dirs=2 -nH "$src" \
				# 		&& pandoc  "${s}.html" -t markdown_github -o README.md \
				# 				&&				 echo OK
				(echo -e "Copyright John Burkardt unless indicated differently.\n\nAll rights reserved\n\nEvery subfolder is published under its own License found in the respective folder's *.html files" \
										 >  LICENSE) \
				&& git add -A															 \
				&& git commit -m "Corrected LICENSE file" \
		    && git gc																		 \
				&& git push origin master
				)
)
done


