>&2 echo -e "\n################################################## stderr for before finish script ##################################################\n"

### bioobject
wget -O bioobject.py https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.1.0/versions.py
python3 bioobject.py

### task_name change "_" to "-"
task_name=$(echo "$task_name" | sed -r 's/_/-/g')

### meta_data
wget -O meta.py https://gitlab.com/intelliseq/workflows/raw/dev/src/main/scripts/bco/1.1.0/meta.py
python3 meta.py "https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/$task_name/$task_version/$task_name.wdl"
jq -r --arg TASK_NAME_WITH_INDEX "$task_name_with_index" '.name = $TASK_NAME_WITH_INDEX' meta.json >meta.json.tmp && mv meta.json.tmp meta.json

### description_domain
task_info=$(jq '. | with_entries(select(.key | contains("input") | not) | select(.key | contains("output") | not) | select(.key | contains("keywords") | not))' meta.json)
inputs=$(jq '. | [.[keys[] | select(contains("input"))]]' meta.json)
outputs=$(jq '. | [.[keys[] | select(contains("output"))]]' meta.json)
keywords=$(jq -s '{ keywords: map(.keywords.keywords[]) }' meta.json)
echo_keywords=`echo $keywords | jq '.keywords'`
task_info_inputs_outputs=$(jo -p task_info="$task_info" input_list="$inputs" output_list="$outputs")
pipeline_steps=$(jo -a "$task_info_inputs_outputs")
description_domain=$(jo -p keywords="$echo_keywords" pipeline_steps="$pipeline_steps")

### bioobject
cpu=$(lscpu | grep '^CPU(s)' | grep -o '[0-9]*')
memory=$(cat /proc/meminfo | grep MemTotal | grep -o '[0-9]*' |  awk '{ print $1/1024/1024 ; exit}')
docker_size=$(cat docker_size)
finishtime=$(date +%s)
starttime=$(cat starttime)
tasktime=$((finishtime-starttime))

os_name=$(cat /etc/os-release | grep -e "^NAME" | grep -o "\".*\"" | sed 's/"//g')
os_version=$(cat /etc/os-release | grep -e "^VERSION_ID" | grep -o "\".*\"" | sed 's/"//g')
jo -p name=$os_name version=$os_version > software_ubuntu.bco
software_prerequisites=$(tools="[]"; for toolfile in $(ls software*.bco); do tools=$(echo $tools | jq ". + [$(cat $toolfile)]") ; done; echo $tools)
if [ "$(ls -A $datasource*.bco)" ]; then
  external_data_endpoints=$(tools="[]"; for toolfile in $(ls datasource*.bco); do tools=$(echo $tools | jq ". + [$(cat $toolfile)]") ; done; echo $tools)
fi

step=$(jo name="$task_name_with_index" version="$task_version")
provenance_domain=$(jo -p \
    steps=$(jo -a $step) \
)

execution_domain=$(jo -p \
  script="https://gitlab.com/intelliseq/workflows/raw/dev/src/main/wdl/tasks/$task_name/$task_version/$task_name.wdl" \
  script_driver="shell" \
  software_prerequisites="$software_prerequisites" \
  external_data_endpoints="$external_data_endpoints" \
  name="$task_name_with_index" \
)

CPU=$(jo param=cpu value=$cpu name="$task_name_with_index")
MEMORY=$(jo param=memory value=$memory name="$task_name_with_index")
TASKTIME=$(jo param=tasktime value=$tasktime name="$task_name_with_index")
DOCKER=$(jo param=docker_image value=$task_docker name="$task_name_with_index")
DOCKER_SIZE=$(jo param=docker_size value=$docker_size name="$task_name_with_index")

parametric_domain=$(jo -a $CPU $MEMORY $TASKTIME $DOCKER $DOCKER_SIZE)

biocomputeobject=$(jo -p \
  bco_spec_version="https://w3id.org/biocompute/1.3.0/" \
  bco_id="https://intelliseq.com/flow/$(uuidgen)" \
  provenance_domain="$provenance_domain" \
  execution_domain="$execution_domain" \
  parametric_domain="$parametric_domain" \
  description_domain="$description_domain" \
)

echo "$biocomputeobject" > bco.json
