from __future__ import division, print_function

import os


def InitialiseWebfolder(webfolder): 
	os.system("mkdir -p "+webfolder)
	webenginesource = "/eos/home-m/mhuwiler/software/php-plots/"
	os.system("cp -r "+webenginesource+"res "+webfolder)
	os.system("cp "+webenginesource+"index.php "+webfolder)
	with open(webenginesource+"example/htaccess", "r") as permissionfile: 
		content = permissionfile.read()
		content = content.replace("/<me>/<my-project>/", webfolder)
		file = open(webfolder+".htaccess", "w")
		file.write(content)
		file.close()
	os.system("cp "+webfolder+".htaccess "+webfolder+"htaccess")


def PublishToWeb(folder, webfolder): 
	webbasepath = "/eos/home-m/mhuwiler/www/Analysis/"
	publicationfolder = "{}{}/".format(webbasepath, webfolder)
	print(publicationfolder)
	InitialiseWebfolder(publicationfolder)
	os.system("cp -r {} {}".format(folder, publicationfolder))
	publicationurl = "https://mhuwiler.web.cern.ch/Analysis/{}/?match=&depth=3".format(webfolder)
	print("Published plots in '{}' to: \n{}".format(folder, publicationurl))

