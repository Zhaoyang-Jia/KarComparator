from jinja2 import Environment, FileSystemLoader
import os

# Define the data
title = "My Text and Images"
text = "This is a sample paragraph with some text."
images = ["/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster0_rotated.png",
          "/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster1_rotated.png",
          "/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster2_rotated.png"]  # List of image file paths

# Create an environment for Jinja2
env = Environment(loader=FileSystemLoader('.'))
template = env.get_template('template.html')

# Render the template with the data
rendered_html = template.render(title=title, text=text, images=images)

# Write the rendered HTML to a file
output_file = 'html_reports/test.html'
with open(output_file, 'w') as f:
    f.write(rendered_html)

print(f"HTML file generated: {os.path.abspath(output_file)}")