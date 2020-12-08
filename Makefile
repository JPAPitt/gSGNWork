add: 
	git add *
	
commit: 
	git commit -m "Update Through MakeFile"

push: 
	git push -u origin master
	
	
all : add commit push
