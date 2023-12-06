import subprocess
import os
import argparse
import re
import signal

# Global variable to keep track of the current downloading file
current_downloading_file = None

def download_files(curl_file, dst):
    global current_downloading_file

    # Check if the download directory exists, and create it if not
    if not os.path.exists(dst):
        os.makedirs(dst)

    # Read the contents of the .curl file
    with open(curl_file, "r") as file:
        curl_commands = file.readlines()

    # Initialize a list to store failed downloads
    failed_downloads = []

    # Loop through each curl command and execute it
    for curl_command in curl_commands:
        try:
            # Extract the file name from the curl command
            file_name = extract_file_name(curl_command)

            # Check if the file already exists in the destination
            if os.path.exists(os.path.join(dst, file_name)):
                print(f"File {file_name} already exists. Skipping download.")
                continue

            # Set the current downloading file
            current_downloading_file = os.path.join(dst, file_name)

            # Use subprocess to run the curl command
            subprocess.run(curl_command.strip(), shell=True, check=True, cwd=dst)

            # Reset the current downloading file after successful download
            current_downloading_file = None
        except subprocess.CalledProcessError as e:
            # Handle the error if the download fails
            failed_downloads.append(curl_command.strip())

    # List the downloads that failed
    if failed_downloads:
        print("The following could not be downloaded:")
        for failed_download in failed_downloads:
            print(failed_download)
    else:
        print("All downloads completed successfully.")

def extract_file_name(curl_command):
    # Use a regular expression to extract the file name from the curl command
    match = re.search(r'-o\s+(\S+)', curl_command)
    return match.group(1) if match else None

def signal_handler(signum, frame):
    # Cleanup function to delete the partial download
    if current_downloading_file and os.path.exists(current_downloading_file):
        print(f"\nScript interrupted. Deleting partial download: {current_downloading_file}")
        os.remove(current_downloading_file)
    exit(1)

def main():
    # Register signal handlers
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Create an argument parser
    parser = argparse.ArgumentParser(description="Download files specified in a .curl file")

    # Add arguments for the .curl file and download directory
    parser.add_argument('--curl_file', required=True, help="Path to the .curl file")
    parser.add_argument('--dst', default=os.getcwd(), help="Path to the download directory (default: current directory)")

    # Parse the command-line arguments
    args = parser.parse_args()

    download_files(args.curl_file, args.dst)

if __name__ == "__main__":
    main()