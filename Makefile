all: pdfs/exercises.pdf pdfs/reading_list.pdf

clean:
	rm -f pdfs/*pdf

pdfs/%.pdf: md/%.md
	pandoc -o $@ $<
