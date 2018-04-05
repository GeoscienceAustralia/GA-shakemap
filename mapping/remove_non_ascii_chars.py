# remove_non_ascii_chars.py

'''
# code tips from:
	https://stackoverflow.com/questions/20078816/replace-non-ascii-characters-with-a-single-space
'''
# load file
text = open('cities1000.txt').read()

# replace chars
newtext = ''.join([i if ord(i) < 128 else ' ' for i in text])

f = open('cities1000_ascii.txt', 'wb')
f.write(newtext)
f.close()

