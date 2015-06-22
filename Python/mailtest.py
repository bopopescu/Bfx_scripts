#!/usr/bin/env python
# encoding: utf-8


def sendMail():
	import smtplib
	from email.mime.text import MIMEText
	fp = open('/home/mevans/logs/hivdb_update_log.txt', 'rb')
	msg = MIMEText(fp.read())	# read logfile into message body
	fp.close()
	me = 'gamera_server@monogrambio.com'
	you = 'bioinformatics@monogrambio.com'
	msg['Subject'] = 'Email test: Do not reply to this address'
	msg['From'] = me
	msg['To'] = you
	# Send the message via our own SMTP server, but don't include the envelope header.
	s = smtplib.SMTP('localhost')
	s.sendmail(me, [you], msg.as_string())
	s.quit()
	return


def main():
	print "\ngenerating email...\n"
	sendMail()
	print "mail sent.\n"

if __name__ == '__main__':
	main()