Étapes pour faire fonctionner en local:


1. Télécharger le projet git

2. Installer une version compilée des drivers sur l'ordinateur en fonction du système d'exploitation en place.
	ALSA (fonctionne sur la majorité des systèmes Linux):
	https://www.alsa-project.org/wiki/Download  -- télécharger alsa-lib, et installer avec `./configure`, puis `make`

3. À l'intérieur du fichier CMake du projet Notable, changer la valeur du champ pour le driver utilisé si le pilote audio diffère des options listées ci-bas:
- ALSA sur Linux
- DirectSound pour Windows
- *Rien à changer pour Mac* (à moins d'exceptions où vous auriez vous-même installé des pilotes)

4. Compiler le projet git de Notable

5. Cela devrait fonctionner. Sinon, il y a probablement eu un problème dans l'installation des pilotes audios à l'étape 2.
