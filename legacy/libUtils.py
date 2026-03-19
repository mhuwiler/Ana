#!/usr/bin/env python



def PromptYesNo(answerasbool=False): 
		# Inspired from Fabrice Couderc 
		rep = ''
		while not rep in [ 'yes', 'no' ]:
			rep = input( "(type 'yes' or 'no'): " ).lower()
		if (answerasbool): 
			if (rep == 'yes'): 
				return True
			else: 
				return False
		return rep


def HoldUntilKeyPress(otherThanEnter=False): 
	# Inspired from Fabrice Couderc 
	rep = ''
	if otherThanEnter: 
		while rep in [ "" ]:
			rep = raw_input( "Press a key to continue... " ).lower()
	else:
		rep = raw_input( "Press a key to continue... " ).lower()
	
	return


def GetDate(invert=True): 
	from datetime import datetime
	date = datetime.today().strftime("%d_%m_%y") 
	if invert: 
		date = date = datetime.today().strftime("%y_%m_%d") 
	return date


def GetFullDate(): 
	from datetime import datetime
	date = datetime.today().strftime("%d_%m_%y__%H_%M_%S")
	if invert: 
		date = datetime.today().strftime("%y_%m_%d__%H_%M_%S")
		return date



