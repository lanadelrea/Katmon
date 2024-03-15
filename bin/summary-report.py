import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

aaf_plot_muts = sys.argv[1]
aaf_plot_amp = sys.argv[2]

# Generate HTML file
html_content = f"""
<!DOCTYPE html>
<html>
<head>
  <title>Alternative Allele Fraction Plot on Mutations Present</title>
</head>
<body>
  <h1>Delta and Omicron specific mutations across the SARS-CoV-2 genome</h1>
  <img src="aaf_plot_muts" alt="aaf per mutation">

  <h1>Alternative Allele Fraction plot per amplicon</h1>
  <img src="aaf_plot_amp" alt="aaf per amplicon">
</body>
</html>
"""

# Write HTML content to a file
with open('plot.html', 'w') as f:
    f.write(html_content)

# Load HTML file and display it
with open('plot.html', 'r') as f:
    html = f.read()
    soup = BeautifulSoup(html, 'html.parser')
    print(soup.prettify())