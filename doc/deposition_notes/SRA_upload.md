Before uploading, organize your files and abide by these recommendations:

    Data upload limits:
        Limit total data uploads to <5 TB total at a time. Your submission subfolder is a temporary holding area and it will be removed once the whole submission is processed.
        Complete a submission through to "Processed" status before uploading more data. If you upload more than the recommended amount of data your files may be removed from our system without warning.
        If you need to submit more than 10 TB data, please contact us first.
    File preparation:
        Compressed fastq files should be <100 GB each.
        Limit your file and folder names to contain only ASCII plain text upper-case and lower-case English letters from A to Z, numbers from 0-9, hyphen, and underscore.
    Directories:
        All files for a submission must be uploaded into a single preload folder. Each preload folder can only be assigned to one SRA submission.
        Do not create complex directory structures of directories within directories.
        Do not upload the same data twice to different directories. If you do this, your data may be removed from our submission system without warning.
        Do not rename a directory after you have linked it to a submission.
        Do not modify or move the files after you have submitted.
    Other recommendations:
        New submitters: Try smaller submissions first to become acquainted with the submission process.
        Please contact us if you have questions.




 Aspera command line provides a fast upload connection over intercontinental distances and has a retry capability for failed connections. Support from your systems administration may be required to set this up properly.
Uploading your data via Aspera command line:

    Download and install the Aspera Connect software (https://www.ibm.com/products/aspera/downloads). From that link click on 'IBM Aspera Connect', then choose the 'Download now...' link for your operating system and install the downloaded app. The ascp program can be found in your install directory: .../Aspera Connect/bin/ascp , C:\Users\\AppData\Local\Programs\IBM\Aspera Connect\bin, etc.
    Download the key file (NOTE: I have saved this here: env\aspera.openssh). Do not include this information or key file on a public page. This information is regularly updated as per our security policies.

    You may use the following command to upload files via Aspera Command-Line:
    ascp -i <path/to/key_file> -QT -l100m -k1 -d <path/to/folder/containing files> subasp@upload.ncbi.nlm.nih.gov:uploads/manfredi.carta_gmail.com_wNN3ofjN

    Where:
    <path/to/key_file> must be an absolute path even if the key is in your current working directory, e.g.: /home/keys/aspera.openssh
    <path/to/folder/containing files> needs to specify the local folder that contains all of the files to upload.
    Using the -d option in your command will ensure this directory is created in your preload area.

    Note: If you upload your files to your account-level directory, you will not be able to use them in your submission. Files uploaded to the account level will be moved to a directory named moved_by_ncbi and subject to deletion after 30 days of inactivity. Therefore, please transfer files to a sub-directory inside your account-level directory.

After upload, select your preload folder from the Files step of your submission.

Note: it takes at least 10 minutes for uploaded files to become available for selection within a submission.

Please complete your submission within 30 days of creating a folder. If you upload files and do not submit them, they will be automatically deleted 30 days after folder creation. 



 You can upload your files via FTP command line or with an FTP client. We recommend using Aspera if you have trouble preloading files with FTP since Aspera can resume interrupted uploads.
Uploading your data via FTP:

    Navigate to the source folder where the files for submission are.
    Establish an FTP connection using the credentials below:
    Address: ftp-private.ncbi.nlm.nih.gov
    Username: subftp
    Password: B1ckyadixdaycsuD
    Please do not share these log-in credentials. Do not include these log-in credentials on a public page. These credentials are changed regularly, as per our security policies.
    SFTP upload is also supported with the same login credentials, but the address should be updated to:
    Address: sftp-private.ncbi.nlm.nih.gov
    Navigate to your account folder:

    From the command line use the 'cd' command:
    cd uploads/manfredi.carta_gmail.com_5SpCDNvD

    When using a GUI-based FTP client (eg: Filezilla, etc.), after you’ve connected to the FTP server, paste your account folder (uploads/manfredi.carta_gmail.com_5SpCDNvD) into the “Remote Site” or “Remote Directory” box on the interface and press “Enter”.

    Until you do this, you will see a message stating "550 /: Permission denied" or "Failed to read the directory listing". We prevent directory listing in the default sign-in folder for security reasons.
    Create a subfolder (required!) with a meaningful name:
    From the command line, use the 'mkdir' (make a directory) command to create a directory/folder with a meaningful name. Here is an example:
    mkdir SUB123456_related_data

    In many GUI-based FTP clients, right-clicking your account folder will give an option to create a sub-directory.

    Note: Use ASCII text only for folder and file names.
    Navigate to the target folder you just created:
    From the command line use the 'cd' command:
    cd SUB123456_related_data

    When using a GUI-based FTP client, click on the new sub-directory to navigate there.
    Copy your files into the target folder:
    From the command line use the 'put' or 'mput' command: put file_name

    In many GUI-based FTP clients, you may drag the local files to the new sub-directory you created in your account folder.

    Note: If you upload your files to your account-level directory, you will not be able to use them in your submission. Files uploaded to the account level will be moved to a directory named moved_by_ncbi and subject to deletion after 30 days of inactivity. Therefore, please transfer files to a sub-directory inside your account-level directory.

After upload, select your preload folder from the Files step of your submission.

Note: it takes at least 10 minutes for uploaded files to become available for selection within a submission.

Please complete your submission within 30 days of creating a folder. If you upload files and do not submit them, they will be automatically deleted 30 days after folder creation. 



More help with using FTP to transfer files: https://www.ncbi.nlm.nih.gov/sra/docs/submitfiles/#file-transfer-protocol-ftp