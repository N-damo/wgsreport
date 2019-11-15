function icon_append(){
    var $=layui.jquery;
    $('h1').prepend('<i class="layui-icon layui-icon-snowflake" style="font-size:30px;margin-right:15px;color:skyblue"></i>');

}



function table_function(){
    var table=layui.table;
    table.init('demo', {height: 500, limit: 10,toolbar:true,page:true,limit:10,defaultToolbar:['filter','exports'],});
    table.render();
}
    
function base_content(){
    var form=layui.form;
    form.render('select');
        $('.base_content-quater').hide();
        var sample=$('#base_content').val();
        $('#base_content-'+sample).show();
        form.on('select(base_content)',function(data){
        $('.base_content-quater').hide();
            var sample=$('#base_content').val();
            $('#base_content-'+sample).show();
        })
    }

function base_quality(){
    var form=layui.form;
    form.render('select');
    $('.base_quality-quater').hide();
    var sample=$('#base_quality').val();
    $('#base_quality-'+sample).show();
    form.on('select(base_quality)',function(data){
    $('.base_quality-quater').hide();
        var sample=$('#base_quality').val();
        $('#base_quality-'+sample).show();
    })
    }


function insert_size(){
    var form=layui.form;
    form.render('select');
    $('.insert_size-quater').hide();
    var sample=$('#insert_size').val();
    $('#insert_size-'+sample).show();
    form.on('select(insert_size)',function(data){
    $('.insert_size-quater').hide();
        var sample=$('#insert_size').val();
        $('#insert_size-'+sample).show();
    })
    }

function depth(){
        var form=layui.form;
        form.render('select');
        $('.depth-quater').hide();
        var sample=$('#depth').val();
        $('#depth-'+sample).show();
        form.on('select(depth)',function(data){
        $('.depth-quater').hide();
            var sample=$('#depth').val();
            $('#depth-'+sample).show();
        })
        }

function genome(){
    var form=layui.form;
    form.render('select');
    $('.genome-quater').hide();
    var sample=$('#genome').val();
    $('#genome-'+sample).show();
    form.on('select(genome)',function(data){
    $('.genome-quater').hide();
        var sample=$('#genome').val();
        $('#genome-'+sample).show();
    })
    }


function shortvariant_stat(){
    var form=layui.form;
    form.render();
    $('#snp').show();
    $('#indel').hide();
    form.on('radio(variant_stat)',function(data){
        if (data.value == 'snp'){
            $('#snp').show();
            $('#indel').hide();
        }
        else {
            $('#snp').hide();
            $('#indel').show();
        }
    })
}

function shortvariant_quality(){
    var form=layui.form;
    form.render();
    $('.snp-quater').hide();
    $('.indel-quater').hide();
    var sample=$('#sample').val();
    $('#snp-'+sample).show();
    form.on('radio()',function(data){
        var variant=data.value;
        var sample=$('#sample').val();
        if (variant == 'snp'){
            $('#snp-'+sample).show();
            $('.indel-quater').hide();
        }
        else {
            $('#indel-'+sample).show();
            $('.snp-quater').hide();
        }
        $('#'+variant+'-'+sample).show();
    })

    form.on('select(sample)',function(data){
        var sample=data.value;
        //alert('重新选择SNP或者INDEL');
        var variant=$('#variant_select input:radio:checked').val();
        $('.snp-quater').hide();
        $('.indel-quater').hide();
        if (variant == 'snp'){
            $('#snp-'+sample).show();
        }
        else {
            $('#indel-'+sample).show();
        }
        $('#'+variant+'-'+sample).show();
        //alert(variant);
        form.on('radio(variant_check)',function(data){
            var variant=data.value;
            $('.snp-quater').hide();
            $('.indel-quater').hide();
            if (variant == 'snp'){
                $('#snp-'+sample).show();
            }
            else {
                $('#indel-'+sample).show();
            }

        })
    })
}

function shortvariant_func(){
    var form=layui.form;
    form.render();
    $('#snp_func').show();
    $('#indel_func').hide();
    form.on('radio(variant_check)',function(data){
        if (data.value == 'snp'){
            $('#snp_func').show();
            $('#indel_func').hide();
        }
        else {
            $('#snp_func').hide();
            $('#indel_func').show();
        }
    })
}


function sv_function(){
    var form=layui.form;
    form.render();
    var sample=$('#sample_selected').val();
    $('.sv-quater').hide();
    $('#sv-'+sample).show();
    form.on('select(cnv)',function(data){
        $('.sv-quater').hide();
        $('#sv-'+data.value).show();
    })
}



function cnv_function(){
    var form=layui.form;
    form.render();
    var sample=$('#sample_selected').val();
    $('.cnv-quater').hide();
    $('#cnv-'+sample).show();
    form.on('select(cnv)',function(data){
        $('.cnv-quater').hide();
        $('#cnv-'+data.value).show();
    })
}

function marker_function(){
    var form=layui.form;
    form.render();
        $('.marker-quater').hide();
        var sample=$('#sample_selected').val();
        $('#marker-'+sample).show();
        form.on('select(sample_selected)',function(data){
        $('.marker-quater').hide();
            var sample=data.value;
            $('#marker-'+sample).show();
        })
    }

function circos_function(){
        var form=layui.form;
        form.render();
            $('.circos-quater').hide();
            var sample=$('#circos').val();
            $('#circos-'+sample).show();
            form.on('select(circos)',function(data){
            $('.circos-quater').hide();
                var sample=data.value;
                $('#circos-'+sample).show();
            })
        }