.PHONY: all build test clean install setupdb

BINARY_DIR=bin
DB_DIR=db
LDFLAGS=-ldflags "-s -w"

all: build

build:
	go build $(LDFLAGS) -o $(BINARY_DIR)/fuse ./cmd/fuse
	go build $(LDFLAGS) -o $(BINARY_DIR)/fuse-getdb ./cmd/fuse-getdb

test:
	go test ./... -v

clean:
	rm -rf $(BINARY_DIR)

install: build
	cp $(BINARY_DIR)/fuse /usr/local/bin/
	cp $(BINARY_DIR)/fuse-getdb /usr/local/bin/

setupdb: build
	$(BINARY_DIR)/fuse --setupdb --datadir $(DB_DIR)