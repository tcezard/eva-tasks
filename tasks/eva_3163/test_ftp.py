import os
from ftplib import FTP
username = 'eva-box-16'
password = '******'
file_to_upload = 'EVA_Submission_template.V1.1.4.xlsx'
destination = 'new_metadata.xlsx'


def upload_file(session, file_to_upload):
    file_name = os.path.basename(file_to_upload)
    with open(file_to_upload, 'rb') as open_file:
        print(session.storbinary(f'STOR {file_name}', open_file))


def list_current_directory(session):
    print(session.retrlines('NLST'))


def download_file(session, file_to_upload, destination):
    file_name = os.path.basename(file_to_upload)
    with open(destination, 'wb') as open_file:
        print(session.retrbinary(f'RETR {file_name}', open_file.write))


ftp = FTP('ftp-private.ebi.ac.uk', timeout=600)
ftp.login(username, password)
# Move to the directory created for upload only
ftp.cwd('test/random_uid')
# Upload works -- 226 Transfer complete.
upload_file(ftp, file_to_upload)
# Listing directory does not work -- 226 Transfer done (but failed to open directory).
list_current_directory(ftp)
# Download same file to new name works -- 226 Transfer complete.
download_file(ftp, file_to_upload, destination)
