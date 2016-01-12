# Download the file

Write-Host ("`nInstalling FastProject to "+(get-location).path+"\FastProject")
Write-Host "`nDownloading Miniconda Installer...`n"

(New-Object System.Net.WebClient).DownloadFile("https://repo.continuum.io/miniconda/Miniconda-latest-Windows-x86_64.exe", "$pwd\Miniconda_Install.exe")

Write-Host "Installing Miniconda...`n"

Start-Process Miniconda_Install.exe "/S /AddToPath=0 /D=$pwd\FastProject" -NoNewWindow -Wait

$env:Path = "$pwd\FastProject\Scripts;" + $env:Path

Write-Host "Installing FastProject dependencies...`n"

conda install numpy scipy matplotlib scikit-learn pandas -y

Write-Host "Installing FastProject...`n"

$fp_archive = (get-childitem FastProject-*.*.*.tar.gz).Name

pip install $fp_archive

# Move fastproject executable to an isolated folder
$fp_script_folder = "$pwd\FastProject\PathScripts"
New-Item $pwd\FastProject\PathScripts -type directory | Out-Null
Move-Item $pwd\FastProject\Scripts\fastproject.exe $fp_script_folder

# Ask user if they want to update path
$title = "Update Path"
$message = "`nDo you want to add the fastproject script to your User PATH?"

$yes = New-Object System.Management.Automation.Host.ChoiceDescription "&Yes", `
    "Prepends the User PATH variable with the location of the FastProject script"

$no = New-Object System.Management.Automation.Host.ChoiceDescription "&No", `
    "User PATH is not modified"

$options = [System.Management.Automation.Host.ChoiceDescription[]]($yes, $no)

$result = $host.ui.PromptForChoice($title, $message, $options, 0) 

if($result -eq 0)
{
    # Update the user's path
    $old_path = (Get-ItemProperty -Path HKCU:\Environment -Name PATH).Path
    $new_path = $fp_script_folder + ";" + $old_path
    cmd /c "setx PATH $new_path"
    Set-ItemProperty -Path HKCU:\Environment -Name PATH -Value $new_path
    Write-Host "User PATH has been updated"
    Write-Host "Restart your terminal to see change"
}
else
{
    Write-Host "User PATH was not modified.`n"
    Write-Host "You may want to add the fastproject script to your path."
    Write-Host "It is located in: $fp_script_folder`n"
}

Write-Host "`nFastProject Successfully Installed"
