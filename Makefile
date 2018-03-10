
all:
	jbuilder build @install

doc: _build/default/_doc/index.html
	jbuilder build @doc

clean:
	jbuilder clean

#	/bin/rm -Rf _build *.install
