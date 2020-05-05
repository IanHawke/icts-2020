all: pdfs/exercises.pdf pdfs/handout.pdf

clean:
	rm -f pdfs/*pdf

pdfs/%.pdf: md/%.md
	pandoc -o $@ $<
