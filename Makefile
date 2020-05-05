all: pdfs/exercises.pdf pdfs/reading_list.pdf

clean:
	rm -f pdfs/*pdf

pdfs/%.pdf: md/%.md
	pandoc --pdf-engine=xelatex -V urlcolor=NavyBlue -V mainfont='Calibri' -V fontsize=12pt -o $@ $<
