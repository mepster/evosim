#/bin/sh

docker exec -it tljh-dev /bin/bash
docker exec --user jupyter-admin -it tljh-dev /bin/bash
docker exec --user jupyter-user1 -it tljh-dev /bin/bash


---


# make users on jh as admin
# login as each user to set pw

docker exec -it tljh-dev /bin/bash
[now on docker]
export JU=jupyter-user7
cd
cp .bashrc .profile /home/$JU
cp -rp /home/jupyter-admin/evosim /home/$JU 
chown -R $JU /home/$JU
chgrp -R $JU /home/$JU


---

# login as root, or a user
#JU=root ; WD=/root ; docker exec --user $JU -w $WD -it tljh-dev /bin/bash
JU=jupyter-user1 ; WD=/home/$JU ; docker exec --user $JU -w $WD -it tljh-dev /bin/bash

# # set up user directories on tljh as root
# JU=root ; WD=/root ; docker exec --user $JU -w $WD -it tljh-dev /bin/bash
# # now on docker
# i=1; while [ "$i" -le 10 ]; do
#     export JU="jupyter-user${i}"
#     echo $JU
#     cp .bashrc .profile /home/$JU
#     cp -rp /home/jupyter-admin/evosim /home/$JU 
#     chown -R $JU /home/$JU
#     chgrp -R $JU /home/$JU
#     usermod -aG conda $JU
#     i=$((i+1))
#     done

# reset all users' evosim repos
# on laptop
export CMD="git pull ; git checkout -f"
for id in admin $(seq 1 10); do
    if [ "$id" = "admin" ]; then
        JU="jupyter-admin"
    else
        JU="jupyter-user${id}"
    fi
    echo "$JU"
    docker exec --user "$JU" -it -w "/home/$JU/evosim" tljh-dev /bin/bash -exec "$CMD"
done

# run any CMD as each user
# on laptop
#export CMD="ls"
#export CMD="/opt/miniforge3/bin/conda env list"
export CMD="git checkout -f ; git pull ; git checkout -f"
for id in admin $(seq 1 10); do
    if [ "$id" = "admin" ]; then
        JU="jupyter-admin"
    else
        JU="jupyter-user${id}"
    fi
    echo "$JU"
    docker exec --user "$JU" -it -w "/home/$JU/evosim" tljh-dev /bin/bash -exec "$CMD"
done
