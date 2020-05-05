all: pdfs/exercises.pdf

clean:
	rm -f pdfs/*pdf

pdfs/%.pdf: md/%.md
	pandoc -o $@ $<
