# SoundEx snippet 1
# from http://code.activestate.com/recipes/52213/
# by Greg Jorgensen 03.06.2001
#
# Warning: This is designed for English names. 
# Warning: This algorithm (by Odell and Russell, as reported in Knuth) is designed for English language surnames. 
# If you have a significant number of non-English surnames, you might do well to alter the values in digits to improve your matches. 
# For example, to accommodate a large number of Spanish surname data, you should count 'J' and 'L' ('L' because of the way 'll' is used) as vowels, setting their position in digit to '0'.
# 
# The basic assumptions of Soundex are that the consonants are more important than the vowels, and that the consonants are grouped into "confusable" groups. 
# Coming up with a set of confusables for a language is not horribly tough, but remember: each group should contain all letters that are confusable with any of those in the group. 
# A slightly better code for both English and Spanish names has digits = '01230120002055012623010202'.
#
# Similar words will return the same soundex
# >>> soundex('hello')
# 'H400'
# >>> soundex('hola')
# 'H400'


def soundex(name, len=4):
    """ soundex module conforming to Knuth's algorithm
        implementation 2000-12-24 by Gregory Jorgensen
        public domain
    """

    # digits holds the soundex values for the alphabet
    digits = '01230120022455012623010202'
    sndx = ''
    fc = ''

    # translate alpha chars in name to soundex digits
    for c in name.upper():
        if c.isalpha():
            if not fc: fc = c   # remember first letter
            d = digits[ord(c)-ord('A')]
            # duplicate consecutive soundex digits are skipped
            if not sndx or (d != sndx[-1]):
                sndx += d

    # replace first digit with first alpha character
    sndx = fc + sndx[1:]

    # remove all 0s from the soundex code
    sndx = sndx.replace('0','')

    # return soundex code padded to len characters
    return (sndx + (len * '0'))[:len]
    


# SoundEx snippet 2
# from http://code.rkevin.com/2010/02/an-implementation-of-the-soundex-algorithm-in-python/
# by Kevin 02.14.2010
#
# Following are the steps to implement the algorithm:
# Ignore all characters in the string being encoded except for the English letters, A to Z.
# The first letter of the Soundex code is the first letter of the string being encoded.
# After the first letter in the string, do not encode vowels or the letters H, W and Y.
# Assign a numeric digit between one and six to all letters except the first using the following mappings:
# 1: B, F, P or V
# 2: C, G, J, K, Q, S, X, Z
# 3: D, T
# 4: L
# 5: M, N
# 6: R
# Where any adjacent digits are the same, remove all but one of those digits unless a vowel, H, W or Y was found between them in the original text.
# Force the code to be four characters in length by padding with zero characters or by truncation
#
# Sample output:
# Smith	S530
# Smythe	S530
# Schultz S243
# Shultz	S432

def get_soundex(name):
	"""Get the soundex code for the string"""
	name = name.upper()

	soundex = ""
	soundex += name[0]

	dictionary = {"BFPV": "1", "CGJKQSXZ":"2", "DT":"3", "L":"4", "MN":"5", "R":"6", "AEIOUHWY":"."}

	for char in name[1:]:
		for key in dictionary.keys():
			if char in key:
				code = dictionary[key]
				if code != soundex[-1]:
					soundex += code

	soundex = soundex.replace(".", "")
	soundex = soundex[:4].ljust(4, "0")

	return soundex

if __name__ == '__main__':
	list = ["Smith", "Smythe", "Robert", "Rupert", "Schultz", "Shultz"]

	print("NAMEttSOUNDEX")
	for name in list:
		print("%stt%s" % (name, get_soundex(name)))
		
