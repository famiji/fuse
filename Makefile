.PHONY: all build test clean install release setupdb

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

# Install the prebuilt local binaries into ~/go/bin (respects GOBIN).
install: build
	@mkdir -p "$${GOBIN:-$$HOME/go/bin}"
	cp $(BINARY_DIR)/fuse $${GOBIN:-$$HOME/go/bin}/fuse
	cp $(BINARY_DIR)/fuse-getdb $${GOBIN:-$$HOME/go/bin}/fuse-getdb

# Cross-compile and upload multi-platform binaries to a GitHub Release.
# Requires GoReleaser (https://goreleaser.com): `brew install goreleaser`.
release:
	goreleaser release --clean

setupdb: build
	$(BINARY_DIR)/fuse --setupdb --datadir $(DB_DIR)