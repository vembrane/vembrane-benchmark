#!env bash
git clone "https://github.com/lindenb/jvarkit.git"
cd jvarkit
git reset --hard 3c11c6de9c11c1ebc5f8c17e19ec304f913c1bee
./gradlew vcffilterjdk
cp dist/vcffilterjdk.jar $CONDA_PREFIX/bin
cat <<EOF > $CONDA_PREFIX/bin/vcffilterjdk
#!/bin/bash
java -jar $CONDA_PREFIX/bin/vcffilterjdk.jar "\$@"
EOF
chmod +x $CONDA_PREFIX/bin/vcffilterjdk